// Copyright 2010, Andreas Biegert

#ifndef CS_PAIRWISE_ALIGNER_INL_H_
#define CS_PAIRWISE_ALIGNER_INL_H_

#include "pairwise_aligner.h"

namespace cs {

template<class Abc>
PairAlignment PairwiseAligner<Abc>::Align(AlignmentMatrices<Abc>& mat) const {
    // Run forward backward
    Forward(mat);
    Backward(mat);

    // Calculate posterior probabilities for MM, MI, IM
    for (size_t i = 0; i <= mat.x.size(); ++i) {
        for (size_t j = 0; j <= mat.y.size(); ++j) {
            mat.P[i][j].MM = mat.F[i][j].MM * mat.B[i][j].MM / mat.ptot;
            mat.P[i][j].MI = mat.F[i][j].MI * mat.B[i][j].MI / mat.ptot;
            mat.P[i][j].IM = mat.F[i][j].IM * mat.B[i][j].IM / mat.ptot;
        }
    }

    // Calculate MAC alignment on posterior matrices
    return MacAlign(mat);
}


template<class Abc>
PairAlignment PairwiseAligner<Abc>::Realign(const PairAlignment& ali, AlignmentMatrices<Abc>& mat) const {
    // Find minimum posterior prob along alignment path and store path in set
    // for quicker access later.
    double pmin = 1.0;
    for (ConstPathIter it = ali.path.begin(); it != ali.path.end(); ++it) {
        if (it->prob < pmin) pmin = it->prob;
    }

    // Subtract 'pmin' from all alignment cells in posterior matrix
    for (ConstPathIter it = ali.path.begin(); it != ali.path.end(); ++it) {
        if (it->pstate == MM) {
            mat.P[it->i][it->j].MM -= pmin;
            if (mat.P[it->i][it->j].MM < 0.0) mat.P[it->i][it->j].MM = 0.0;
        }
        // switch (it->pstate) {
        //     case MM:
        //         mat.P[it->i][it->j].MM -= pmin;
        //         if (mat.P[it->i][it->j].MM < 0.0) mat.P[it->i][it->j].MM = 0.0;
        //         break;
        //     case MI:
        //         mat.P[it->i][it->j].MI -= pmin;
        //         if (mat.P[it->i][it->j].MI < 0.0) mat.P[it->i][it->j].MI = 0.0;
        //         break;
        //     case IM:
        //         mat.P[it->i][it->j].IM -= pmin;
        //         if (mat.P[it->i][it->j].IM < 0.0) mat.P[it->i][it->j].IM = 0.0;
        //         break;
        // }
    }

    // Calculate MAC alignment on altered posterior matrix
    PairAlignment newali(MacAlign(mat));
    newali.subopt = ali.subopt + 1;

    return newali;
}


// Calculates co-emission probability between abstract-state columns 'px' and 'py'
// normalized with abstract-state background frequencies in array 'bg'
template<class Abc>
inline double PairwiseAligner<Abc>::ContextProb(const SparseProfileCol& px, const SparseProfileCol& py) const {
    if (context_score_ == 0.0 || (px.empty() && py.empty())) {
        return 1.0;
    } else if (px.num_nonempty() < py.num_nonempty()) {
        double sc = 0.0;
        for (SparseProfileCol::ElementIter xk = px.begin(); xk != px.end(); ++xk) {
            sc += (xk->prob * py[xk->index] * cons_[xk->index]) / priors_[xk->index];
        }
        LOG(DEBUG) << "contextprob=" << pow(sc + 1.0, context_score_);
        return pow(sc + 1.0, context_score_);
    } else {
        double sc = 0.0;
        for (SparseProfileCol::ElementIter yk = py.begin(); yk != py.end(); ++yk) {
            sc += (yk->prob * px[yk->index] * cons_[yk->index]) / priors_[yk->index];
        }
        LOG(DEBUG) << "contextprob=" << pow(sc + 1.0, context_score_);
        return pow(sc + 1.0, context_score_);
    }
    return 1.0;
}

// template<class Abc>
// inline double PairwiseAligner<Abc>::ContextProb(const SparseProfileCol& px, const SparseProfileCol& py) const {
//     if (context_score_ == 0.0 || (px.empty() && py.empty())) {
//         return 1.0;
//     } else if (px.num_nonempty() < py.num_nonempty()) {
//         double sc = 0.0, xi = 0.0, yj = 0.0;
//         for (SparseProfileCol::ElementIter xk = px.begin(); xk != px.end(); ++xk) {
//             sc += (xk->prob * py[xk->index]) / priors_[xk->index];
//             xi += xk->prob;
//         }
//         for (SparseProfileCol::ElementIter yk = py.begin(); yk != py.end(); ++yk)
//             yj += yk->prob;
//         LOG(DEBUG) << "contextprob=" << pow(sc + 1.0 - xi * yj, context_score_);
//         return pow(sc + 1.0 - xi * yj, context_score_);
//     } else {
//         double sc = 0.0, xi = 0.0, yj = 0.0;
//         for (SparseProfileCol::ElementIter yk = py.begin(); yk != py.end(); ++yk) {
//             sc += (yk->prob * px[yk->index]) / priors_[yk->index];
//             yj += yk->prob;
//         }
//         for (SparseProfileCol::ElementIter xk = px.begin(); xk != px.end(); ++xk)
//             xi += xk->prob;
//         LOG(DEBUG) << "contextprob=" << pow(sc + 1.0 - xi * yj, context_score_);
//         return pow(sc + 1.0 - xi * yj, context_score_);
//     }
//     return 1.0;
// }

// Calculates phylogeny match probability between profile columns 'px' and 'py'
template<class Abc>
inline double PairwiseAligner<Abc>::MatchProb(const ProfileColumn<Abc>& px, const ProfileColumn<Abc>& py) const {
    double rv = 0.0;

    if (maty_ == NULL) {  // column score
        for (size_t a = 0; a < Abc::kSize; ++a)
            rv += (px[a] * py[a]) / matx_->p(a);

    } else {  // phylogeny aware match score
        static double bg[] = { 0.3, 0.2, 0.2, 0.3 };  // FIXME: these background freqs work only for DNA!
        // Calculate model probability
        for (size_t a = 0; a < Abc::kSize; ++a) {
            double sum_px = 0.0;
            for (size_t b = 0; b < Abc::kSize; ++b)
                sum_px += matx_->r(b,a) * px[b];
            double sum_py = 0.0;
            for (size_t b = 0; b < Abc::kSize; ++b)
                sum_py += maty_->r(b,a) * py[b];
            rv += bg[a] * sum_px * sum_py;
        }
        // Calculate null-model probability
        double null_px = 0.0;
        for (size_t a = 0; a < Abc::kSize; ++a)
            null_px += bg[a] * px[a];
        double null_py = 0.0;
        for (size_t a = 0; a < Abc::kSize; ++a)
            null_py += bg[a] * py[a];
        rv /= null_px * null_py;
    }

    return rv;
}


template<class Abc>
void PairwiseAligner<Abc>::Forward(AlignmentMatrices<Abc>& mat) const {
    const Graph& x = mat.x.g;
    const Graph& y = mat.y.g;
    const size_t xL = mat.x.size();
    const size_t yL = mat.y.size();
    const VertexVec& x_vertices = mat.x.vertices;
    const VertexVec& y_vertices = mat.y.vertices;

    ProbMatrix& F = mat.F;
    Vector<double>& scale = mat.scale;
    InEdgeIter ei, ej, ei_end, ej_end;
    OutEdgeIter ek, ek_end;
    VertexIndex index;
    size_t i, j, ii, jj;
    double pmax_i;

    // Initialization of top left cell
    F[0][0].MM = 1.0; F[0][0].MI = F[0][0].IM = 0.0;

    // Initialization of top row, i.e. cells (0,j)
    scale[0] = 1.0;
    for (size_t l = 1; l <= yL; ++l) {
        j = y_vertices[l];
        F[0][j].MM = F[0][j].MI = F[0][j].IM = 0.0;
        for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
            jj = index[source(*ej, y)];
            F[0][j].IM += y[*ej].weight *
                ( F[0][jj].MM * x[0].M2I * y[jj].M2M +   // MM -> IM
                  F[0][jj].IM * x[0].I2I * y[jj].M2M );  // IM -> IM
        }
    }

    // Core dynamic programming
    for (size_t k = 1; k <= xL; ++k) {
        i = x_vertices[k];
        pmax_i = 0;

        // Initialize leftmost cell (i,0)
        F[i][0].MM = F[i][0].MI = F[i][0].IM = 0.0;
        for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
            ii = index[source(*ei, x)];
            F[i][0].MI += scale[i] * x[*ei].weight *
                ( F[ii][0].MM * x[ii].M2M * y[0].M2I +   // MM -> MI
                  F[ii][0].MI * x[ii].M2M * y[0].I2I );  // MI -> MI
        }

        // Fill remaining cells in row i
        for (size_t l = 1; l <= yL; ++l) {
            j = y_vertices[l];

            // MM recursion relation
            F[i][j].MM = 0.0;
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                    jj = index[source(*ej, y)];

                    F[i][j].MM += scale[i] * x[*ei].weight * y[*ej].weight *
                        ( F[ii][jj].MM * x[ii].M2M * y[jj].M2M +   // MM -> MM
                          F[ii][jj].MI * x[ii].M2M * y[jj].I2M +   // MI -> MM
                          F[ii][jj].IM * x[ii].I2M * y[jj].M2M );  // IM -> MM
                }
            }
            F[i][j].MM *= MatchProb(x[i].probs, y[j].probs) * ContextProb(x[i].contexts, y[j].contexts);

            // MI recursion relation
            F[i][j].MI = 0.0;
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                F[i][j].MI += scale[i] * x[*ei].weight *
                    ( F[ii][j].MM * x[ii].M2M * y[j].M2I +   // MM -> MI
                      F[ii][j].MI * x[ii].M2M * y[j].I2I );  // MI -> MI
            }

            // IM recursion relation
            F[i][j].IM = 0.0;
            for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                jj = index[source(*ej, y)];
                F[i][j].IM += y[*ej].weight *
                    ( F[i][jj].MM * x[i].M2I * y[jj].M2M +   // MM -> IM
                      F[i][jj].IM * x[i].I2I * y[jj].M2M );  // IM -> IM
            }

            if (F[i][j].MM > pmax_i) pmax_i = F[i][j].MM;
        }

        // Calculate scale factors for target vertices of vertex i
        for (tie(ek, ek_end) = out_edges(i, x); ek != ek_end; ++ek) {
            ii = index[target(*ek, x)];
            if (ii != kStartEndVertex && 1.0 / (pmax_i + 1.0) < scale[ii])
                scale[ii] = 1.0 / (pmax_i + 1.0);
        }
    }

    // Calculate total probability by summing over all transitions into END-state
    double ptot = 0.0;
    for (tie(ei, ei_end) = in_edges(kStartEndVertex, x); ei != ei_end; ++ei) {
        i = index[source(*ei, x)];
        for (tie(ej, ej_end) = in_edges(kStartEndVertex, y); ej != ej_end; ++ej) {
            j = index[source(*ej, y)];
            ptot += F[i][j].MM * x[i].M2M * y[j].M2M * x[*ei].weight * y[*ej].weight;
            ptot += F[i][j].MI * x[i].M2M * y[j].I2M * x[*ei].weight * y[*ej].weight;
            ptot += F[i][j].IM * x[i].I2M * y[j].M2M * x[*ei].weight * y[*ej].weight;
        }
    }
    mat.ptot = ptot;
    LOG(INFO) << "ptot=" << ptot;
}


template<class Abc>
void PairwiseAligner<Abc>::Backward(AlignmentMatrices<Abc>& mat) const {
    const Graph& x = mat.x.g;
    const Graph& y = mat.y.g;
    const size_t xL = mat.x.size();
    const size_t yL = mat.y.size();
    const VertexVec& x_vertices = mat.x.vertices;
    const VertexVec& y_vertices = mat.y.vertices;

    ProbMatrix& B = mat.B;
    Vector<double>& scale = mat.scale;
    OutEdgeIter ei, ej, ei_end, ej_end;
    InEdgeIter ek, el, ek_end, el_end;
    VertexIndex index;
    size_t i, j, ii, jj;

    // Recursion relations
    for (int k = xL; k >= 0; --k) {
        i = x_vertices[k];
        for (int l = yL; l >= 0; --l) {
            j = y_vertices[l];
            B[i][j].MM = B[i][j].MI = B[i][j].IM = 0.0;

            for (tie(ei, ei_end) = out_edges(i, x); ei != ei_end; ++ei) {
                ii = index[target(*ei, x)];
                for (tie(ej, ej_end) = out_edges(j, y); ej != ej_end; ++ej) {
                    jj = index[target(*ej, y)];
                    double pmatch = x[*ei].weight * y[*ej].weight * scale[ii];

                    if (ii != kStartEndVertex && jj != kStartEndVertex)
                        pmatch *= B[ii][jj].MM * MatchProb(x[ii].probs, y[jj].probs) * ContextProb(x[ii].contexts, y[jj].contexts);
                    else if (ii != kStartEndVertex || jj != kStartEndVertex)
                        continue;

                    B[i][j].MM += pmatch * x[i].M2M * y[j].M2M;
                    B[i][j].MI += pmatch * x[i].M2M * y[j].I2M;
                    B[i][j].IM += pmatch * x[i].I2M * y[j].M2M;
                }
            }

            for (tie(ei, ei_end) = out_edges(i, x); ei != ei_end; ++ei) {
                ii = index[target(*ei, x)];
                if (ii == kStartEndVertex) continue;
                double tmp = x[*ei].weight * B[ii][j].MI * scale[ii];
                B[i][j].MM += tmp * x[i].M2M * y[j].M2I;
                B[i][j].MI += tmp * x[i].M2M * y[j].I2I;
            }

            for (tie(ej, ej_end) = out_edges(j, y); ej != ej_end; ++ej) {
                jj = index[target(*ej, y)];
                if (jj == kStartEndVertex) continue;
                double tmp = y[*ej].weight * B[i][jj].IM;
                B[i][j].MM += tmp * x[i].M2I * y[j].M2M;
                B[i][j].IM += tmp * x[i].I2I * y[j].M2M;
            }
        }
    }
}


template<class Abc>
PairAlignment PairwiseAligner<Abc>::MacAlignNormalized(AlignmentMatrices<Abc>& mat) const {
    const Graph& x = mat.x.g;
    const Graph& y = mat.y.g;
    const size_t xL = mat.x.size();
    const size_t yL = mat.y.size();
    const VertexVec& x_vertices = mat.x.vertices;
    const VertexVec& y_vertices = mat.y.vertices;

    ProbMatrix& P = mat.P;
    ScoreMatrix& S = mat.S;
    BacktraceMatrix& b = mat.b;
    InEdgeIter ei, ej, ei_end, ej_end;
    VertexIndex index;
    size_t i, j, ii, jj, ii_max, jj_max;
    uint8_t XX_max;
    ScoreCell sc_max, sc;

    // Initialization of top row, i.e. cells (0,j)
    S[0][0] = ScoreCell(); b[0][0] = BacktraceCell(0, 0, MM);
    for (size_t l = 1; l <= yL; ++l) {
        j = y_vertices[l];
        jj_max = 0; sc_max.avg = 0.0;
        for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
            jj = index[source(*ej, y)];
            sc.sum = S[0][jj].sum + gapf_ * P[0][j].IM;
            sc.len = S[0][jj].len + 1;
            sc.avg = sc.sum / sc.len;
            if (sc.avg > sc_max.avg) { sc_max = sc; jj_max = jj; }
        }
        S[0][j] = sc_max;
        b[0][j] = BacktraceCell(0, jj_max, IM);
    }

    // Initialization of leftmost column, i.e. cells (i,0)
    for (size_t k = 1; k <= xL; ++k) {
        i = x_vertices[k];
        ii_max = 0; sc_max.avg = 0.0;
        for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
            ii = index[source(*ei, x)];
            sc.sum = S[ii][0].sum + gapf_ * P[i][0].MI;
            sc.len = S[ii][0].len + 1;
            sc.avg = sc.sum / sc.len;
            if (sc.avg > sc_max.avg) { sc_max = sc; ii_max = ii; }
        }
        S[i][0] = sc_max;
        b[i][0] = BacktraceCell(ii_max, 0, MI);
    }

    // Core dynamic programming
    for (size_t k = 1; k <= xL; ++k) {
        i = x_vertices[k];
        for (size_t l = 1; l <= yL; ++l) {
            j = y_vertices[l];
            ii_max = 0; jj_max = 0; XX_max = MM; sc_max.avg = 0.0;
            // Maximimize over all MM possibilities
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                    jj = index[source(*ej, y)];
                    sc.sum = S[ii][jj].sum + P[i][j].MM;
                    sc.len = S[ii][jj].len + 1;
                    sc.avg = sc.sum / sc.len;
                    if (sc.avg > sc_max.avg) {
                        sc_max = sc;
                        ii_max = ii;
                        jj_max = jj;
                        XX_max = MM;
                    }
                }
            }
            // Maximize over all MI possibilities
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                sc.sum = S[ii][j].sum + gapf_ * P[i][j].MI;
                sc.len = S[ii][j].len + 1;
                sc.avg = sc.sum / sc.len;
                if (sc.avg > sc_max.avg) {
                    sc_max = sc;
                    ii_max = ii;
                    jj_max = j;
                    XX_max = MI;
                }
            }
            // Maximimize over all IM possibilities
            for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                jj = index[source(*ej, y)];
                sc.sum = S[i][jj].sum + gapf_ * P[i][j].IM;
                sc.len = S[i][jj].len + 1;
                sc.avg = sc.sum / sc.len;
                if (sc.avg > sc_max.avg) {
                    sc_max = sc;
                    ii_max = i;
                    jj_max = jj;
                    XX_max = IM;
                }
            }

            // Set S[i][j] to calculated maximum and store backtrace pointer
            S[i][j] = sc_max;
            b[i][j] = BacktraceCell(ii_max, jj_max, XX_max);
        }
    }

    // Find maximum score of  "corner-like" cells, where both vertices transit
    // into END-state. From there we need to start our backtrace.
    ii_max = 0; jj_max = 0;
    double avg_max = 0.0;
    for (tie(ei, ei_end) = in_edges(kStartEndVertex, x); ei != ei_end; ++ei) {
        ii = index[source(*ei, x)];
        for (tie(ej, ej_end) = in_edges(kStartEndVertex, y); ej != ej_end; ++ej) {
            jj = index[source(*ej, y)];
            if (S[ii][jj].avg > avg_max) {
                avg_max = S[ii][jj].avg;
                ii_max = ii;
                jj_max = jj;
            }
        }
    }

    // Run backtrace starting from 'ii_max', 'jj_max'
    i = ii_max; j = jj_max;
    PairAlignment ali;
    while (i != 0 || j != 0) {
        double prob = 0.0;
        switch (b[i][j].pstate) {
            case MM:
                prob = P[i][j].MM;
                ali.num_matches++;
                break;
            case MI:
                prob = P[i][j].MI;
                break;
            case IM:
                prob = P[i][j].IM;
                break;
        }

        // Add cell to alignment path
        ali.path.push_back(Step(i, j, b[i][j].pstate, prob));
        ali.sum_probs += prob;
        ali.pmin = MIN(ali.pmin, prob);

        // Take care to when advancing to predecessor cell
        ii = b[i][j].from_i;
        jj = b[i][j].from_j;

        i = ii;
        j = jj;
    }
    ali.avg_probs = ali.sum_probs / ali.path.size();
    reverse(ali.path.begin(), ali.path.end());

    return ali;
}


template<class Abc>
PairAlignment PairwiseAligner<Abc>::MacAlign(AlignmentMatrices<Abc>& mat) const {
    const Graph& x = mat.x.g;
    const Graph& y = mat.y.g;
    const size_t xL = mat.x.size();
    const size_t yL = mat.y.size();
    const VertexVec& x_vertices = mat.x.vertices;
    const VertexVec& y_vertices = mat.y.vertices;

    ProbMatrix& P = mat.P;
    ScoreMatrix& S = mat.S;
    BacktraceMatrix& b = mat.b;
    InEdgeIter ei, ej, ei_end, ej_end;
    VertexIndex index;
    size_t i, j, ii, jj, ii_max, jj_max;
    uint8_t XX_max;
    ScoreCell sc_max, sc;

    // Initialization of top row, i.e. cells (0,j)
    S[0][0] = ScoreCell(); b[0][0] = BacktraceCell(0, 0, MM);
    for (size_t l = 1; l <= yL; ++l) {
        j = y_vertices[l];
        jj_max = 0; sc_max.sum = 0.0;
        for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
            jj = index[source(*ej, y)];
            sc.sum = S[0][jj].sum + gapf_ * P[0][j].IM * y[j].col.size();
            if (sc.sum > sc_max.sum) { sc_max = sc; jj_max = jj; }
        }
        S[0][j] = sc_max;
        b[0][j] = BacktraceCell(0, jj_max, IM);
    }

    // Initialization of leftmost column, i.e. cells (i,0)
    for (size_t k = 1; k <= xL; ++k) {
        i = x_vertices[k];
        ii_max = 0; sc_max.sum = 0.0;
        for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
            ii = index[source(*ei, x)];
            sc.sum = S[ii][0].sum + gapf_ * P[i][0].MI * x[i].col.size();
            if (sc.sum > sc_max.sum) { sc_max = sc; ii_max = ii; }
        }
        S[i][0] = sc_max;
        b[i][0] = BacktraceCell(ii_max, 0, MI);
    }

    // Core dynamic programming
    for (size_t k = 1; k <= xL; ++k) {
        i = x_vertices[k];
        for (size_t l = 1; l <= yL; ++l) {
            j = y_vertices[l];
            ii_max = 0; jj_max = 0; XX_max = MM; sc_max.sum = 0.0;
            // Maximimize over all MM possibilities
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                    jj = index[source(*ej, y)];
                    sc.sum = S[ii][jj].sum + P[i][j].MM * (x[i].col.size() + y[j].col.size());
                    if (sc.sum > sc_max.sum) {
                        sc_max = sc;
                        ii_max = ii;
                        jj_max = jj;
                        XX_max = MM;
                    }
                }
            }
            // Maximize over all MI possibilities
            for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
                ii = index[source(*ei, x)];
                sc.sum = S[ii][j].sum + gapf_ * P[i][j].MI * x[i].col.size();
                if (sc.sum > sc_max.sum) {
                    sc_max = sc;
                    ii_max = ii;
                    jj_max = j;
                    XX_max = MI;
                }
            }
            // Maximimize over all IM possibilities
            for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
                jj = index[source(*ej, y)];
                sc.sum = S[i][jj].sum + gapf_ * P[i][j].IM * y[j].col.size();
                if (sc.sum > sc_max.sum) {
                    sc_max = sc;
                    ii_max = i;
                    jj_max = jj;
                    XX_max = IM;
                }
            }

            // Set S[i][j] to calculated maximum and store backtrace pointer
            S[i][j] = sc_max;
            b[i][j] = BacktraceCell(ii_max, jj_max, XX_max);
        }
    }

    // Find maximum score of  "corner-like" cells, where both vertices transit
    // into END-state. From there we need to start our backtrace.
    ii_max = 0; jj_max = 0;
    double sum_max = 0.0;
    for (tie(ei, ei_end) = in_edges(kStartEndVertex, x); ei != ei_end; ++ei) {
        ii = index[source(*ei, x)];
        for (tie(ej, ej_end) = in_edges(kStartEndVertex, y); ej != ej_end; ++ej) {
            jj = index[source(*ej, y)];
            if (S[ii][jj].sum > sum_max) {
                sum_max = S[ii][jj].sum;
                ii_max = ii;
                jj_max = jj;
            }
        }
    }

    // Run backtrace starting from 'ii_max', 'jj_max'
    i = ii_max; j = jj_max;
    PairAlignment ali;
    while (i != 0 || j != 0) {
        double prob = 0.0;
        switch (b[i][j].pstate) {
            case MM:
                prob = P[i][j].MM;
                ali.num_matches++;
                break;
            case MI:
                prob = P[i][j].MI;
                break;
            case IM:
                prob = P[i][j].IM;
                break;
        }

        // Add cell to alignment path
        ali.path.push_back(Step(i, j, b[i][j].pstate, prob));
        ali.sum_probs += prob;
        // ali.pmin = MIN(ali.pmin, prob);
        if (b[i][j].pstate == MM) ali.pmin = MIN(ali.pmin, prob);

        // Take care to when advancing to predecessor cell
        ii = b[i][j].from_i;
        jj = b[i][j].from_j;

        i = ii;
        j = jj;
    }
    ali.avg_probs = ali.sum_probs / ali.path.size();
    reverse(ali.path.begin(), ali.path.end());

    return ali;
}


inline std::string PairStateToString(uint8_t xx) {
    std::string rv;
    switch (xx) {
        case MM: rv = "MM"; break;
        case MI: rv = "MI"; break;
        case IM: rv = "IM"; break;
        default: throw Exception("Unknown pair state '%d'!", xx);
    }
    return rv;
}

inline std::ostream& operator<< (std::ostream& out, const PairAlignment& ali) {
    out << "Alignment " << ali.subopt << std::endl;
    out << "sum-probs:   " << ali.sum_probs << std::endl;
    out << "avg-probs:   " << ali.avg_probs << std::endl;
    out << "pmin:        " << ali.pmin << std::endl;
    out << "num-matches: " << ali.num_matches << std::endl;
    out << "length:      " << ali.path.size() << std::endl;

    out << strprintf("%2s  %4s  %4s  %4s\n", "XX", "i", "j", "Prob");
    for (ConstPathIter it = ali.path.begin(); it != ali.path.end(); ++it) {
        std::string s = PairStateToString(it->pstate);
        out << strprintf("%2s  %4zu  %4zu  %4.2f\n", s.c_str(), it->i, it->j, it->prob);
    }
    return out;
}

template<class Abc>
std::ostream& operator<< (std::ostream& out, const AlignmentMatrices<Abc>& m) {
    // out << "Forward matrix F[i][j]:" << std::endl;

    // out << strprintf("%3s %-2s | ", "", "");
    // for (size_t j = 0; j < m.F.ncols(); ++j)
    //   out << strprintf("%6zu  ", j);
    // out << std::endl << std::string(9 + m.F.ncols() * 8, '-') << std::endl;

    // for (size_t i = 0; i < m.F.nrows(); ++i) {
    //   out << strprintf("%3s %-2s | ", "", "MM");
    //   for (size_t j = 0; j < m.F.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.F[i][j].MM);

    //   out << strprintf("\n%3zu %-2s | ", i, "MI");
    //   for (size_t j = 0; j < m.F.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.F[i][j].MI);

    //   out << strprintf("\n%3s %-2s | ", "", "IM");
    //   for (size_t j = 0; j < m.F.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.F[i][j].IM);

    //   out << std::endl;
    //   if (i + 1 != m.F.nrows()) out << strprintf("%3s %-2s | \n", "", "");
    // }

    // out << "Backward matrix B[i][j]:" << std::endl;

    // out << strprintf("%3s %-2s | ", "", "");
    // for (size_t j = 0; j < m.B.ncols(); ++j)
    //   out << strprintf("%6zu  ", j);
    // out << std::endl << std::string(9 + m.B.ncols() * 8, '-') << std::endl;

    // for (size_t i = 0; i < m.B.nrows(); ++i) {
    //   out << strprintf("%3s %-2s | ", "", "MM");
    //   for (size_t j = 0; j < m.B.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.B[i][j].MM);

    //   out << strprintf("\n%3zu %-2s | ", i, "MI");
    //   for (size_t j = 0; j < m.B.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.B[i][j].MI);

    //   out << strprintf("\n%3s %-2s | ", "", "IM");
    //   for (size_t j = 0; j < m.B.ncols(); ++j)
    //     out << strprintf("%6.1g  ", m.B[i][j].IM);

    //   out << std::endl;
    //   if (i + 1 != m.B.nrows()) out << strprintf("%3s %-2s | \n", "", "");
    // }

    out << "Scaling factors for forward/backward probabilities:" << std::endl;
    for (size_t i = 0; i < m.S.nrows(); ++i)
        out << strprintf("scale[%zu] = %6.4f\n", i, m.scale[i]);

    out << "Posterior probability matrix P[i][j]:" << std::endl;

    out << strprintf("%3s %-2s | ", "", "");
    for (size_t j = 0; j < m.P.ncols(); ++j)
        out << strprintf("%6zu  ", j);
    out << std::endl << std::string(9 + m.P.ncols() * 8, '-') << std::endl;

    for (size_t i = 0; i < m.P.nrows(); ++i) {
        out << strprintf("%3s %-2s | ", "", "MM");
        for (size_t j = 0; j < m.P.ncols(); ++j)
            out << strprintf("%6.2f  ", m.P[i][j].MM * 100);

        out << strprintf("\n%3zu %-2s | ", i, "MI");
        for (size_t j = 0; j < m.P.ncols(); ++j)
            out << strprintf("%6.2f  ", m.P[i][j].MI * 100);

        out << strprintf("\n%3s %-2s | ", "", "IM");
        for (size_t j = 0; j < m.P.ncols(); ++j)
            out << strprintf("%6.2f  ", m.P[i][j].IM * 100);
        out << std::endl;
        if (i + 1 != m.P.nrows()) out << strprintf("%3s %-2s | \n", "", "");
    }

    // out << "Posterior match probability matrix P[i][j].MM:" << std::endl;
    // out << strprintf("%3s | ", "");
    // for (size_t j = 0; j < m.P.ncols(); ++j)
    //   out << strprintf("%3zu ", j);
    // out << std::endl << std::string(6 + (m.P.ncols() - 1) * 4, '-') << std::endl;

    // for (size_t i = 0; i < m.P.nrows(); ++i) {
    //   out << strprintf("%3zu | ", i);
    //   for (size_t j = 0; j < m.P.ncols(); ++j)
    //     out << strprintf("%3d ", iround(m.P[i][j].MM * 100));
    //   out << std::endl;
    // }

    // out << "MAC score matrix S[i][j]" << std::endl;
    // out << strprintf("%5s | ", "");
    // for (size_t j = 0; j < m.S.ncols(); ++j)
    //   out << strprintf("%5zu ", j);
    // out << std::endl << std::string(7 + m.S.ncols() * 6, '-') << std::endl;

    // for (size_t i = 0; i < m.S.nrows(); ++i) {
    //   out << strprintf("%5zu | ", i);
    //   for (size_t j = 0; j < m.S.ncols(); ++j)
    //     out << strprintf("%5.2f ", m.S[i][j].avg);
    //   out << std::endl;
    // }

    // out << "Backtrace matrix b[i][j]" << std::endl;
    // out << strprintf("%5s | ", "");
    // for (size_t j = 0; j < m.b.ncols(); ++j)
    //   out << strprintf("%5zu ", j);
    // out << std::endl << std::string(7 + m.b.ncols() * 6, '-') << std::endl;

    // for (size_t i = 0; i < m.b.nrows(); ++i) {
    //   out << strprintf("%5zu | ", i);
    //   for (size_t j = 0; j < m.b.ncols(); ++j) {
    //     std::string s;
    //     if (m.b[i][j].pstate == MM) s = "MM";
    //     if (m.b[i][j].pstate == MI) s = "^^";
    //     if (m.b[i][j].pstate == IM) s = "<<";
    //     out << strprintf("%5s ", s.c_str());
    //   }
    //   out << std::endl;
    // }

    return out;
}

}  // namespace cs

#endif  // CS_PAIRWISE_ALIGNER_INL_H_
