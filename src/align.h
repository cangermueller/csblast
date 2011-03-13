// Copyright 2009, Andreas Biegert

#ifndef CS_ALIGN_H_
#define CS_ALIGN_H_

#include "po_hmm-inl.h"

namespace cs {

// Pair states needed for backtracing
enum PairState { MM=0, MI=1, IM=2 };

// Many of these make up an alignment path.
struct Path {
  int i, j;
  uint8_t pstate;
};

// Keeps track of path through dynamic programing matrix
struct BacktraceCell {
  uint8_t MM, MI, IM;
};

// Keeps track of forward/backward probs
struct ProbCell {
  ProbCell() : MM(0.0), MI(0.0), IM(0.0) {}
  double MM, MI, IM;
};

// POD with forward backward matrices
struct DynProgMatrices {
  DynProgMatrices(size_t nq, size_t np)
      : F(nq, np), B(nq, np), P(nq, np), ptot(0.0) {}

  friend std::ostream& operator<< (std::ostream& out, const DynProgMatrices& m) {
    out << "Forward matrix F[i][j]:" << std::endl;

    out << strprintf("%3s %-2s | ", "", "");
    for (size_t j = 0; j < m.F.ncols(); ++j)
      out << strprintf("%6zu  ", j);
    out << std::endl << std::string(9 + m.F.ncols() * 8, '-') << std::endl;

    for (size_t i = 0; i < m.F.nrows(); ++i) {
      out << strprintf("%3s %-2s | ", "", "MM");
      for (size_t j = 0; j < m.F.ncols(); ++j)
        out << strprintf("%6.1g  ", m.F[i][j].MM);

      out << strprintf("\n%3zu %-2s | ", i, "MI");
      for (size_t j = 0; j < m.F.ncols(); ++j)
        out << strprintf("%6.1g  ", m.F[i][j].MI);

      out << strprintf("\n%3s %-2s | ", "", "IM");
      for (size_t j = 0; j < m.F.ncols(); ++j)
        out << strprintf("%6.1g  ", m.F[i][j].IM);

      out << std::endl;
      if (i + 1 != m.F.nrows()) out << strprintf("%3s %-2s | \n", "", "");
    }

    out << "Backward matrix B[i][j]:" << std::endl;

    out << strprintf("%3s %-2s | ", "", "");
    for (size_t j = 0; j < m.B.ncols(); ++j)
      out << strprintf("%6zu  ", j);
    out << std::endl << std::string(9 + m.B.ncols() * 8, '-') << std::endl;

    for (size_t i = 0; i < m.B.nrows(); ++i) {
      out << strprintf("%3s %-2s | ", "", "MM");
      for (size_t j = 0; j < m.B.ncols(); ++j)
        out << strprintf("%6.1g  ", m.B[i][j].MM);

      out << strprintf("\n%3zu %-2s | ", i, "MI");
      for (size_t j = 0; j < m.B.ncols(); ++j)
        out << strprintf("%6.1g  ", m.B[i][j].MI);

      out << strprintf("\n%3s %-2s | ", "", "IM");
      for (size_t j = 0; j < m.B.ncols(); ++j)
        out << strprintf("%6.1g  ", m.B[i][j].IM);

      out << std::endl;
      if (i + 1 != m.B.nrows()) out << strprintf("%3s %-2s | \n", "", "");
    }

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

    out << "Posterior match probability matrix P[i][j].MM:" << std::endl;

    out << strprintf("%3s | ", "");
    for (size_t j = 1; j < m.P.ncols(); ++j)
      out << strprintf("%3zu ", j);
    out << std::endl << std::string(6 + (m.P.ncols() - 1) * 4, '-') << std::endl;

    for (size_t i = 1; i < m.P.nrows(); ++i) {
      out << strprintf("%3zu | ", i);
      for (size_t j = 1; j < m.P.ncols(); ++j)
        out << strprintf("%3d ", iround(m.P[i][j].MM * 100));
      out << std::endl;
    }
    return out;
  }

  Matrix<ProbCell> F;
  Matrix<ProbCell> B;
  Matrix<ProbCell> P;
  double ptot;
};

// Calculates co-emission probability between abstract-state columns 'px' and 'py'
// normalized with abstract-state background frequencies in array 'bg'
inline double CoEmission(const SparseProfileCol& px,
                             const SparseProfileCol& py,
                             const double* bg) {
  double rv = 0.0;
  for (SparseProfileCol::ElementIter xk = px.begin(); xk != px.end(); ++xk)
    rv += xk->prob * py[xk->index] / bg[xk->index];
  return rv;
}

template<class Abc>
double Forward(const POHmm<Abc>& hmm_x,
               const POHmm<Abc>& hmm_y,
               const double* bg,
               Matrix<ProbCell>& F) {
  typedef typename POHmm<Abc>::Graph Graph;
  typedef typename POHmm<Abc>::InEdgeIter InEdgeIter;
  typedef typename POHmm<Abc>::VertexIndex VertexIndex;
  typedef typename POHmm<Abc>::VertexVector IndexVector;

  const Graph& x = hmm_x.g;
  const Graph& y = hmm_y.g;
  const size_t xL = hmm_x.size();
  const size_t yL = hmm_y.size();
  IndexVector x_indices(hmm_x.GetSortedVertices());
  IndexVector y_indices(hmm_y.GetSortedVertices());
  InEdgeIter ei, ej, ei_end, ej_end;
  VertexIndex index;
  size_t i, j, ii, jj;

  assert_eq(0, static_cast<int>(x_indices[0]));
  assert_eq(0, static_cast<int>(y_indices[0]));

  // Initialization of top row, i.e. cells (0,j)
  F[0][0].MM = 1.0; F[0][0].MI = F[0][0].IM = 0.0;
  for (size_t l = 1; l <= yL; ++l) {
    j = y_indices[l];
    F[0][j].MM = F[0][j].MI = F[0][j].IM = 0.0;
    for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
      jj = index[source(*ej, y)];
      F[0][j].IM += y[*ej].prob *
        ( F[0][jj].MM * x[0].tr[M2I] * y[jj].tr[M2M] +   // MM -> IM
          F[0][jj].IM * x[0].tr[I2I] * y[jj].tr[M2M] );  // IM -> IM
    }
  }

  // Initialization of leftmost column, i.e. cells (i,0)
  for (size_t k = 1; k <= xL; ++k) {
    i = x_indices[k];
    F[i][0].MM = F[i][0].MI = F[i][0].IM = 0.0;
    for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
      ii = index[source(*ei, x)];
      F[i][0].MI += x[*ei].prob *
        ( F[ii][0].MM * x[ii].tr[M2M] * y[0].tr[M2I] +   // MM -> MI
          F[ii][0].MI * x[ii].tr[M2M] * y[0].tr[I2I] );  // MI -> MI
    }
  }

  // Core dynamic programming
  for (size_t k = 1; k <= xL; ++k) {
    i = x_indices[k];
    for (size_t l = 1; l <= yL; ++l) {
      j = y_indices[l];

      // MM recursion relation
      F[i][j].MM = 0.0;
      for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
        ii = index[source(*ei, x)];
        for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
          jj = index[source(*ej, y)];

          F[i][j].MM += x[*ei].prob * y[*ej].prob *
            ( F[ii][jj].MM * x[ii].tr[M2M] * y[jj].tr[M2M] +   // MM -> MM
              F[ii][jj].MI * x[ii].tr[M2M] * y[jj].tr[I2M] +   // MI -> MM
              F[ii][jj].IM * x[ii].tr[I2M] * y[jj].tr[M2M] );  // IM -> MM
        }
      }
      F[i][j].MM *= CoEmission(x[i].contexts, y[j].contexts, bg);

      // MI recursion relation
      F[i][j].MI = 0.0;
      for (tie(ei, ei_end) = in_edges(i, x); ei != ei_end; ++ei) {
        ii = index[source(*ei, x)];
        F[i][j].MI += x[*ei].prob *
          ( F[ii][j].MM * x[ii].tr[M2M] * y[j].tr[M2I] +   // MM -> MI
            F[ii][j].MI * x[ii].tr[M2M] * y[j].tr[I2I] );  // MI -> MI
      }

      // IM recursion relation
      F[i][j].IM = 0.0;
      for (tie(ej, ej_end) = in_edges(j, y); ej != ej_end; ++ej) {
        jj = index[source(*ej, y)];
        F[i][j].IM += y[*ej].prob *
          ( F[i][jj].MM * x[i].tr[M2I] * y[jj].tr[M2M] +   // MM -> IM
            F[i][jj].IM * x[i].tr[I2I] * y[jj].tr[M2M] );  // IM -> IM
      }
    }
  }

  // Calculate total probability by summing over all transitions into END-state
  double ptot = 0.0;
  for (tie(ei, ei_end) = in_edges(kStartEndVertex, x); ei != ei_end; ++ei) {
    i = index[source(*ei, x)];
    for (tie(ej, ej_end) = in_edges(kStartEndVertex, y); ej != ej_end; ++ej) {
      j = index[source(*ej, y)];
      ptot += F[i][j].MM * x[i].tr[M2M] * y[j].tr[M2M] * x[*ei].prob * y[*ej].prob;
      ptot += F[i][j].MI * x[i].tr[M2M] * y[j].tr[I2M] * x[*ei].prob * y[*ej].prob;
      ptot += F[i][j].IM * x[i].tr[I2M] * y[j].tr[M2M] * x[*ei].prob * y[*ej].prob;
    }
  }

  LOG(ERROR) << "ptot=" << ptot;
  return ptot;
}


template<class Abc>
void Backward(const POHmm<Abc>& hmm_x,
              const POHmm<Abc>& hmm_y,
              const double* bg,
              Matrix<ProbCell>& B) {
  typedef typename POHmm<Abc>::Graph Graph;
  typedef typename POHmm<Abc>::OutEdgeIter OutEdgeIter;
  typedef typename POHmm<Abc>::InEdgeIter InEdgeIter;
  typedef typename POHmm<Abc>::VertexIndex VertexIndex;
  typedef typename POHmm<Abc>::VertexVector IndexVector;

  const Graph& x = hmm_x.g;
  const Graph& y = hmm_y.g;
  const size_t xL = hmm_x.size();
  const size_t yL = hmm_y.size();
  IndexVector x_indices(hmm_x.GetSortedVertices());
  IndexVector y_indices(hmm_y.GetSortedVertices());
  OutEdgeIter ei, ej, ei_end, ej_end;
  InEdgeIter ek, el, ek_end, el_end;
  VertexIndex index;
  size_t i, j, ii, jj;
  Bitset is_end_x(hmm_x.size() + 1);  // keeps track of before END vertices in x
  Bitset is_end_y(hmm_y.size() + 1);  // keeps track of before END vertices in y

  assert_eq(0, static_cast<int>(x_indices[0]));
  assert_eq(0, static_cast<int>(y_indices[0]));

  // Initialize "corner-like" cells, where both vertices transit into END-state
  for (tie(ek, ek_end) = in_edges(kStartEndVertex, x); ek != ek_end; ++ek) {
    i = index[source(*ek, x)];
    is_end_x.set(i);
    for (tie(el, el_end) = in_edges(kStartEndVertex, y); el != el_end; ++el) {
      j = index[source(*el, y)];
      // FIXME: take edge probs into account, also for total Forward prob!
      B[i][j].MM = x[i].tr[M2M] * y[j].tr[M2M] * x[*ek].prob * y[*el].prob;
      B[i][j].MI = x[i].tr[M2M] * y[j].tr[I2M] * x[*ek].prob * y[*el].prob;
      B[i][j].IM = x[i].tr[I2M] * y[j].tr[M2M] * x[*ek].prob * y[*el].prob;
      is_end_y.set(j);
    }
  }

  // Initialization of bottom row, i.e. cells (xL,j)
  for (tie(ek, ek_end) = in_edges(kStartEndVertex, x); ek != ek_end; ++ek) {
    i = index[source(*ek, x)];
    for (size_t l = yL; l > 0; --l) {
      j = y_indices[l];
      if (!is_end_y.test(j)) {
        B[i][j].MM = B[i][j].MI = B[i][j].IM = 0.0;
        for (tie(ej, ej_end) = out_edges(j, y); ej != ej_end; ++ej) {
          jj = index[target(*ej, y)];
          B[i][j].MM += y[*ej].prob * B[i][jj].IM * x[i].tr[M2I] * y[j].tr[M2M];
          B[i][j].IM += y[*ej].prob * B[i][jj].IM * x[i].tr[I2I] * y[j].tr[M2M];
        }
      }
    }
  }

  // Initialize rightmost column, i.e. cells (i,yL)
  for (tie(el, el_end) = in_edges(kStartEndVertex, y); el != el_end; ++el) {
    j = index[source(*el, y)];
    for (size_t k = xL; k > 0; --k) {
      i = x_indices[k];
      if (!is_end_x.test(i)) {
        B[i][j].MM = B[i][j].MI = B[i][j].IM = 0.0;
        for (tie(ei, ei_end) = out_edges(i, x); ei != ei_end; ++ei) {
          ii = index[target(*ei, x)];
          B[i][j].MM += x[*ei].prob * B[ii][j].MI * x[i].tr[M2M] * y[j].tr[M2I];
          B[i][j].MI += x[*ei].prob * B[ii][j].MI * x[i].tr[M2M] * y[j].tr[I2I];
        }
      }
    }
  }

  // Core dynamic programming
  for (size_t k = xL - 1; k > 0; --k) {
    i = x_indices[k];
    for (size_t l = yL - 1; l > 0; --l) {
      j = y_indices[l];
      B[i][j].MM = B[i][j].MI = B[i][j].IM = 0.0;

      for (tie(ei, ei_end) = out_edges(i, x); ei != ei_end; ++ei) {
        ii = index[target(*ei, x)];
        for (tie(ej, ej_end) = out_edges(j, y); ej != ej_end; ++ej) {
          jj = index[target(*ej, y)];
          double pmatch = x[*ei].prob * y[*ej].prob * B[ii][jj].MM *
            CoEmission(x[ii].contexts, y[jj].contexts, bg);

          B[i][j].MM += pmatch * x[i].tr[M2M] * y[j].tr[M2M];
          B[i][j].MI += pmatch * x[i].tr[M2M] * y[j].tr[I2M];
          B[i][j].IM += pmatch * x[i].tr[I2M] * y[j].tr[M2M];
        }
      }

      for (tie(ei, ei_end) = out_edges(i, x); ei != ei_end; ++ei) {
        ii = index[target(*ei, x)];
        double tmp = x[*ei].prob * B[ii][j].MI;

        B[i][j].MM += tmp * x[i].tr[M2M] * y[j].tr[M2I];
        B[i][j].MI += tmp * x[i].tr[M2M] * y[j].tr[I2I];
      }

      for (tie(ej, ej_end) = out_edges(j, y); ej != ej_end; ++ej) {
        jj = index[target(*ej, y)];
        double tmp = y[*ej].prob * B[i][jj].IM;

        B[i][j].MM += tmp * x[i].tr[M2I] * y[j].tr[M2M];
        B[i][j].IM += tmp * x[i].tr[I2I] * y[j].tr[M2M];
      }
    }
  }
}

}  // namespace cs

#endif  // CS_ALIGN_H_
