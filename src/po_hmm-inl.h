// Copyright 2009, Andreas Biegert

#ifndef CS_PO_HMM_INL_H_
#define CS_PO_HMM_INL_H_

#include "po_hmm.h"
#include "context_library-inl.h"
#include "emission.h"
#include "matrix_pseudocounts-inl.h"
#include "tamura_nei_matrix.h"

namespace cs {

template<class Abc>
POHmm<Abc>::POHmm(double gapo, double gape)
        : num_alis(0), neff(0.0), gap_open(gapo), gap_extend(gape) {
    add_vertex(VertexProperties<Abc>(), g);
}

template<class Abc>
POHmm<Abc>::POHmm(const Sequence<Abc>& seq, double gapo, double gape)
        : seqw(1, 1.0),
          aliw(1, 1.0),
          num_alis(1),
          neff(1.0),
          gap_open(gapo),
          gap_extend(gape) {
    // Add sequence
    seqs.push_back(seq);

    // Add start vertex
    Vertex v = add_vertex(VertexProperties<Abc>(), g);
    g[v].col.push_back(std::make_pair(0, -1));
    g[v].alis.set(0);
    g[v].flow = 1.0;

    // Add one vertex per residue and connect vertices left to right
    size_t nv = 1;
    for (size_t i = 0; i < seq.length(); ++i, ++nv) {
        v = add_vertex(VertexProperties<Abc>(), g);
        g[v].col.push_back(std::make_pair(0, i));
        g[v].alis.set(0);
        g[v].flow = 1.0;
        add_edge(nv - 1, nv, EdgeProperties(1.0,  1.0), g);
    }
    // Connect last vertex to end-vertex
    add_edge(nv - 1, kStartEndVertex, EdgeProperties(1.0,  1.0), g);

    // Initialize remaining members
    Init();
}

template<class Abc>
POHmm<Abc>::POHmm(const Alignment<Abc>& ali, double gapo, double gape)
        : seqw(ali.nseqs(), 0.0),
          aliw(1, 1.0),
          num_alis(1),
          neff(0.0),
          gap_open(gapo),
          gap_extend(gape) {
    const size_t nseqs = ali.nseqs();
    assert_eq(0, ali.ninsert());
    // Add alignment sequences
    for (size_t n = 0; n < nseqs; ++n) seqs.push_back(ali.GetSequence(n));

    // Add start vertex
    Vertex v = add_vertex(VertexProperties<Abc>(), g);
    g[v].alis.set(0);
    g[v].flow = 1.0;
    for (size_t k = 0; k < ali.nseqs(); ++k)
        g[v].col.push_back(std::make_pair(k, -1));

    // Add one vertex per match column and connect vertices left to right as well
    // as adding indel edges for jumping over gaps.
    Vector<size_t> last_res(nseqs, kStartEndVertex);
    size_t nv = 1;
    for (size_t i = 0; i < ali.nmatch(); ++i, ++nv) {
        // Add vertex
        v = add_vertex(VertexProperties<Abc>(), g);
        g[v].alis.set(0);
        g[v].flow = 1.0;
        // Add incoming edges
        Bitset in_edges(ali.nmatch() + 1);
        for (size_t n = 0; n < ali.nseqs(); ++n) {
            if (ali[i][n] < Abc::kGap) {
                if (!in_edges.test(last_res[n])) {
                    add_edge(last_res[n], nv, EdgeProperties(1.0,  1.0), g);
                    in_edges.set(last_res[n]);
                }
                last_res[n] = nv;
            }
        }
    }

    // Connect vertices representing end residues to end-vertex
    Bitset in_edges(ali.nmatch() + 1);
    for (size_t n = 0; n < ali.nseqs(); ++n) {
        if (!in_edges.test(last_res[n])) {
            add_edge(last_res[n], kStartEndVertex, EdgeProperties(1.0,  1.0), g);
            in_edges.set(last_res[n]);
        }
    }

    // Add sequence-residue pairs to col vectors
    for (size_t k = 0; k < ali.nseqs(); ++k) {
        size_t r = 0;
        for (size_t i = 0; i < ali.nmatch(); ++i)
            if (ali.seq(k, i) < Abc::kGap)
                g[i+1].col.push_back(std::make_pair(k, r++));
    }

    // Initialize remaining members
    Init();
}

template<class Abc>
void POHmm<Abc>::AssignCountsAndProbs() {
    // Columns with sum of seq-weights below 'insert_thresh' get Neff=1
    const double insert_thresh = 1.0 / sqrt(2.0 * seqs.size());
    LOG(ERROR) << "insert_thresh=" << insert_thresh;

    // Reset counts in profile columns to zero
    for (VertexVecIter it = vertices.begin() + 1; it != vertices.end(); ++it) {
        Assign(g[*it].counts, 0.0);
        Assign(g[*it].probs, 0.0);
        g[*it].neff = 0.0;
    }

    // Calculate average Neff over all vertices and assign probs and counts
    neff = 0.0;
    double sum_flow = 0.0;
    double f[Abc::kSize];  // residue frequencies in column i
    for (size_t i = 1; i <= size(); ++i) {
        // Reset frequency array
        for (size_t a = 0; a < Abc::kSize; ++a) f[a] = 0.0;
        // Calculate frequencies
        double wi = 0.0;  // sum of seq-weight in column i
        for (ColIter r = g[i].col.begin(); r != g[i].col.end(); ++r) {
            if (seqs[r->first][r->second] < Abc::kAny)
                f[seqs[r->first][r->second]] += seqw[r->first];
            wi += seqw[r->first];
        }
        Normalize(&f[0], Abc::kSize);
        // Calculate Neff contribution of column i
        // if (wi >= insert_thresh) { // this will make all columns in a two sequence alignment matches
        if (wi > insert_thresh) {
            double ni = seqs.size() + 1.0;
            for (size_t a = 0; a < Abc::kSize; ++a)
                ni -= SQR(f[a]) * seqs.size();
            neff += g[i].flow * ni;
            sum_flow += g[i].flow;
        }
        // Assign residue frequencies to probs and counts vector
        for (size_t a = 0; a < Abc::kSize; ++a) {
            g[i].probs[a]  = f[a];
            g[i].counts[a] = f[a];  // counts still need to be normalized!
        }
        g[i].probs[Abc::kAny] = 1.0;
    }
    neff /= sum_flow;
    LOG(ERROR) << "Average Neff=" << neff;

    // Normalize counts either to one if it's an insert column or to Neff otherwise
    for (size_t i = 1; i <= size(); ++i) {
        double wi = 0.0;  // sum of seq-weight in column i
        for (ColIter r = g[i].col.begin(); r != g[i].col.end(); ++r)
            wi += seqw[r->first];
        g[i].neff = wi < insert_thresh ? 1.0 : neff;
        Normalize(g[i].counts, g[i].neff);
    }
}

template<class Abc>
void POHmm<Abc>::AssignTransitionProbs() {
    VertexIter vi, vertex_end;
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        g[*vi].M2M = 1.0 - gap_open;
        g[*vi].M2I = gap_open;
        g[*vi].I2M = 1.0 - gap_extend;
        g[*vi].I2I = gap_extend;
    }
}

// TODO: use dynamic bitset instead of static one so we can get rid of
// kMaxSeqCapacity and hold as many seqs as we like.
template<class Abc>
void POHmm<Abc>::AssignEdgeWeights() {
    typedef std::bitset<kMaxSeqCapacity> SeqBitset;
    Vector<int> v2idx(vertices.size());
    VertexIter vi, vertex_end;

    // Build reverse mapping from vertex to its index in sorted vertex vector
    for (size_t i = 0; i < vertices.size(); ++i)
        v2idx[vertices[i]] = i;

    // Forward pass through out-edges
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        if (out_degree(*vi, g) == 1) {
            // Fast path: there is only one out-edge, so we set its prob to one
            g[*out_edges(*vi, g).first].weight = 1.0;
        } else {
            // Slow path: more than one out edge, we need to calc the edge prob based
            // on global sequence weights and vertex flows
            std::vector<OutEdgeIter> iters;
            OutEdgeIter ei, edge_end;
            for (tie(ei, edge_end) = out_edges(*vi, g); ei != edge_end; ++ei)
                iters.push_back(ei);
            std::sort(iters.begin(), iters.end(), OutEdgeIterCompare(g, v2idx));

            // Calculate sum of sequence weights at source vertex 'v' and default bitmask
            double vw = 0.0;
            SeqBitset vseqs;
            for (ColIter r = g[*vi].col.begin(); r != g[*vi].col.end(); ++r) {
                vw += seqw[r->first];
                vseqs.set(r->first);
            }
            std::vector<SeqBitset> branch_seqs(num_alis, vseqs);

            for (size_t i = 0; i < iters.size(); ++i) {
                Vertex u = target(*iters[i], g);
                // Build seqs bitmask for this vertex 'u'
                SeqBitset useqs(vseqs);
                for (size_t n = 0; n < num_alis; ++n)
                    if (g[u].alis.test(n)) useqs &= branch_seqs[n];
                // Calculate sum of sequence weights at target vertex 'u'
                double uw = 0.0;
                for (ConstColIter r = g[u].col.begin(); r != g[u].col.end(); ++r) {
                    if (useqs.test(r->first)) {
                        uw += seqw[r->first];
                        useqs.set(r->first, false);  // turn sequence off
                    }
                }
                // Calculate branch weight based on flow
                AliBitset vu_alis(g[*vi].alis & g[u].alis);
                double vu_flow = 0.0;
                for (size_t n = 0; n < num_alis; ++n)
                    if (vu_alis.test(n)) vu_flow += aliw[n];
                // Set edge weight
                g[*iters[i]].weight = (vu_flow / g[*vi].flow) * (uw / vw);
                // Update bitmasks of all branches that used vertex 'u'
                for (size_t n = 0; n < num_alis; ++n)
                    if (g[u].alis.test(n)) branch_seqs[n] &= useqs;
            }
        }
    }

    // Backward pass through in-edges
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        if (in_degree(*vi, g) == 1) {
            // Fast path: there is only one out-edge, so we set its prob to one
            g[*in_edges(*vi, g).first].weight_rev = 1.0;
        } else {
            // Slow path: more than one out edge, we need to calc the edge prob based
            // on global sequence weights and vertex flows
            std::vector<InEdgeIter> iters;
            InEdgeIter ei, edge_end;
            for (tie(ei, edge_end) = in_edges(*vi, g); ei != edge_end; ++ei)
                iters.push_back(ei);
            std::sort(iters.begin(), iters.end(), InEdgeIterCompare(g, v2idx));

            // Calculate sum of sequence weights at source vertex 'v' and default bitmask
            double vw = 0.0;
            SeqBitset vseqs;
            for (ColIter r = g[*vi].col.begin(); r != g[*vi].col.end(); ++r) {
                vw += seqw[r->first];
                vseqs.set(r->first);
            }
            std::vector<SeqBitset> branch_seqs(num_alis, vseqs);

            for (size_t i = 0; i < iters.size(); ++i) {
                Vertex u = source(*iters[i], g);
                // Build seqs bitmask for this vertex 'u'
                SeqBitset useqs(vseqs);
                for (size_t n = 0; n < num_alis; ++n)
                    if (g[u].alis.test(n)) useqs &= branch_seqs[n];
                // Calculate sum of sequence weights at source vertex 'u'
                double uw = 0.0;
                for (ConstColIter r = g[u].col.begin(); r != g[u].col.end(); ++r) {
                    if (useqs.test(r->first)) {
                        uw += seqw[r->first];
                        useqs.set(r->first, false);  // turn sequence off
                    }
                }
                // Calculate branch weight based on flow
                AliBitset vu_alis(g[*vi].alis & g[u].alis);
                double vu_flow = 0.0;
                for (size_t n = 0; n < num_alis; ++n)
                    if (vu_alis.test(n)) vu_flow += aliw[n];
                // Set edge weight
                g[*iters[i]].weight_rev = (vu_flow / g[*vi].flow) * (uw / vw);
                // Update bitmasks of all branches that used vertex 'u'
                for (size_t n = 0; n < num_alis; ++n)
                    if (g[u].alis.test(n)) branch_seqs[n] &= useqs;
            }
        }
    }
}

template<class Abc>
void POHmm<Abc>::AssignAncestralProbs(const POHmm<Abc>& x,
                                      const POHmm<Abc>& y,
                                      const SubstitutionMatrix<Abc>& matx,
                                      const SubstitutionMatrix<Abc>& maty) {
    VertexIter vi, vertex_end;
    VertexIndex index;
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex) {
            if (g[*vi].pstate == MI) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    g[*vi].probs[a] = x.g[g[*vi].i].probs[a];
            } else if (g[*vi].pstate == IM) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    g[*vi].probs[a] = y.g[g[*vi].j].probs[a];
            } else {
                for (size_t a = 0; a < Abc::kSize; ++a) {
                    double sum_px = 0.0;
                    for (size_t b = 0; b < Abc::kSize; ++b)
                        sum_px += matx.r(b,a) * x.g[g[*vi].i].probs[b];
                    double sum_py = 0.0;
                    for (size_t b = 0; b < Abc::kSize; ++b)
                        sum_py += maty.r(b,a) * y.g[g[*vi].j].probs[b];
                    g[*vi].probs[a] = sum_px * sum_py;
                }
            }
        }
    }
}

template<class Abc>
void POHmm<Abc>::DisableSkipped(size_t max_skip_count) {
    VertexIter vi, vertex_end;
    VertexIndex index;
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex) {
            if (g[*vi].skip_count > max_skip_count) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    g[*vi].probs[a] = 0.0;
            }
        }
    }
}

template<class Abc>
void POHmm<Abc>::AssignContextStateProbs(const ContextLibrary<Abc>& lib, const Emission<Abc>& emission, double cutoff) {
    VertexIter vi, vertex_end;
    VertexIndex index;
    Vector<double> pp(lib.size());
    for (tie(vi, vertex_end) = boost::vertices(g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex) {
            CalculatePosteriorProbs(lib, emission, g, *vi, &pp[0]);
            g[*vi].contexts = SparseProfileCol(&pp[0], lib.size(), cutoff);
        }
    }
}

template<class Abc>
void POHmm<Abc>::CalculateSequenceWeights() {
    // Calculate global seq-weights
    seqw.assign(seqs.size(), 0.0);
    size_t c[Abc::kSizeAny];  // number of times a residue occurs at column i
    for (size_t i = 1; i <= size(); ++i) {
        // Reset residue counts
        for (size_t a = 0; a < Abc::kSizeAny; ++a) c[a] = 0;
        // Count residues
        for (ColIter r = g[i].col.begin(); r != g[i].col.end(); ++r)
            c[seqs[r->first][r->second]]++;
        // Count number of  different residues
        size_t na = 0;
        for (size_t a = 0; a < Abc::kSize; ++a) if (c[a]) na++;
        if (na == 0) na = 1;  // col consists of only gaps and ANYs
        // Add seq-weight contribution
        for (ColIter r = g[i].col.begin(); r != g[i].col.end(); ++r)
            if (seqs[r->first][r->second] < Abc::kAny)
                seqw[r->first] += g[i].flow / (c[seqs[r->first][r->second]] * na * seqs[r->first].length());
    }
    Normalize(&seqw[0], seqs.size());
}

template<class Abc>
void POHmm<Abc>::Init() {
    // Topologically sort PO-HMM vertices
    typedef filtered_graph<Graph, EdgeFilter> FilteredGraph;
    EdgeFilter filter(&g);
    FilteredGraph fg(g, filter);
    vertices.clear();
    topological_sort(fg, std::back_inserter(vertices));
    reverse(vertices.begin(), vertices.end());
    assert_eq(kStartEndVertex, static_cast<size_t>(vertices.front()));

    // Calculate index of each vertex in topolical sorted vertex vector
    indices.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) indices[vertices[i]] = i;

    // Init remaining members and property structs
    CalculateSequenceWeights();
    AssignEdgeWeights();
    AssignTransitionProbs();
    AssignCountsAndProbs();
}

template<class Abc>
Alignment<Abc> POHmm<Abc>::GetAlignment(size_t n) const {
    // Determine number of alignment columns
    size_t ncols = 0;
    for (ConstVertexVecIter it = vertices.begin() + 1; it != vertices.end(); ++it)
        if (g[*it].alis.test(n)) ncols++;

    // Construct the alignment by filling in aligned residues at each column
    Alignment<Abc> ali(ncols, seqs.size());
    size_t i = 0;
    for (ConstVertexVecIter it = vertices.begin() + 1; it != vertices.end(); ++it) {
        if (g[*it].alis.test(n)) {
            for (ConstColIter r = g[*it].col.begin(); r != g[*it].col.end(); ++r)
                ali[i][r->first] = seqs[r->first][r->second];
            i++;  // advance column index
        }
    }

    // Set sequence headers
    for (size_t k = 0; k < seqs.size(); ++k)
        ali.set_header(k, seqs[k].header());

    return ali;
}

template<class Abc>
void POHmm<Abc>::SanityCheck() const {
    for (size_t n = 0; n < num_alis; ++n) {
        Alignment<Abc> ali(GetAlignment(n));
        for (size_t k = 0; k < seqs.size(); ++k) {
            Sequence<Abc> ali_seq(ali.GetSequence(k));
            if (seqs[k].length() != ali_seq.length()) {
                LOG(ERROR) << "PO-HMM alignment:\n" << ali;
                LOG(ERROR) << "Original sequence:\n" << seqs[k];
                LOG(ERROR) << "PO-HMM sequence:\n" << ali_seq;
                throw Exception("PO-HMM sanity check: sequence '%s' in alignment %zu should have %zu residues but actually has %zu",
                                ali_seq.header().c_str(), n, seqs[k].length(), ali_seq.length());
            }
            for (size_t i = 0; i < ali_seq.length(); ++i) {
                if (seqs[k][i] != ali_seq[i]) {
                    LOG(ERROR) << "PO-HMM alignment:\n" << ali;
                    LOG(ERROR) << "Original sequence:\n" << seqs[k];
                    LOG(ERROR) << "PO-HMM sequence:\n" << ali_seq;
                    throw Exception("PO-HMM sanity check: sequence '%s' in alignment %zu differs at position %zu", ali_seq.header().c_str(), n, i);
                }
            }
        }
    }
}

// Prints a PO-HMM in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const POHmm<Abc>& hmm) {
    typedef typename POHmm<Abc>::ConstVertexVecIter ConstVertexVecIter;
    typedef typename POHmm<Abc>::Graph Graph;

    const Graph& g = hmm.g;
    const size_t slen = 5;
    const size_t sidx = kMaxAliCapacity - slen;
    typename POHmm<Abc>::VertexIndex index = get(vertex_index, hmm.g);
    typename POHmm<Abc>::OutEdgeIter ei, edge_end;

    // Construct multiple alignment on the fly for printing of ali-columns with gaps
    Alignment<Dna> ali(hmm.size(), hmm.seqs.size());
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it) {
        for (ConstColIter r = g[*it].col.begin(); r != g[*it].col.end(); ++r) {
            ali[*it - 1][r->first] = hmm.seqs[r->first][r->second];
        }
    }

    out << "PO-HMM:" << std::endl;
    // Print vertex indices
    out << strprintf("%-10s", "Vertex");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4d  ", static_cast<int>(index[*it]));
    out << std::endl << std::endl;

    // Print alignment
    for (size_t n = 0; n < hmm.seqs.size(); ++n) {
        out << strprintf("%-10s", hmm.seqs[n].header().substr(0, 10).c_str());
        for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
            out << strprintf("   %c  ", ali.chr(n, *it - 1));
        out << std::endl;
    }
    out << std::endl;

    // Print alignment bits
    out << strprintf("%-10s", "Alis");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%5s ", g[*it].alis.to_string().substr(sidx, slen).c_str());
    out << std::endl << std::endl;

    // Print flow
    out << strprintf("%-10s", "Flow");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].flow);
    out << std::endl << std::endl;

    // Print neff
    out << strprintf("%-10s", "Neff");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].neff);
    out << std::endl << std::endl;

    // Print profile counts
    for (size_t a = 0; a < Abc::kSize; ++a) {
        out << strprintf("Counts %-3c", Abc::kIntToChar[a]);
        for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
            out << strprintf("%4.2f  ", g[*it].counts[a]);
        out << std::endl;
    }
    out << std::endl;

    // Print profile probs
    for (size_t a = 0; a < Abc::kSize; ++a) {
        out << strprintf("Probs %-4c", Abc::kIntToChar[a]);
        for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
            out << strprintf("%4.2f  ", g[*it].probs[a]);
        out << std::endl;
    }
    out << std::endl;

    // Print transition probs
    out << strprintf("Trans %-4s", "M2M");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].M2M);
    out << strprintf("\nTrans %-4s", "M2I");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].M2I);
    out << strprintf("\nTrans %-4s", "I2M");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].I2M);
    out << strprintf("\nTrans %-4s", "I2I");
    for (ConstVertexVecIter it = hmm.begin() + 1; it != hmm.end(); ++it)
        out << strprintf("%4.2f  ", g[*it].I2I);
    out << std::endl << std::endl;

    // Print edge-weights
    out << "Edges:" << std::endl;
    for (ConstVertexVecIter it = hmm.begin(); it != hmm.end(); ++it) {
        out << strprintf("%3d (%4.2f, %zu, %4s)  ",
                         index[*it], g[*it].flow, g[*it].col.size(),
                         g[*it].alis.to_string().substr(sidx, slen).c_str());
        for (tie(ei, edge_end) = out_edges(*it, g); ei != edge_end; ++ei) {
            out << strprintf("   == %4.2f, %4.2f ==> %-3d", g[*ei].weight,
                             g[*ei].weight_rev, index[target(*ei, g)]);
        }
        out << std::endl;
    }
    out << std::endl;

    out << "Alignment weights:" << std::endl;
    out << StringifyRange(&hmm.aliw[0], &hmm.aliw[0] + hmm.aliw.size());
    out << std::endl;
    out << "Global sequence weights:" << std::endl;
    out << StringifyRange(&hmm.seqw[0], &hmm.seqw[0] + hmm.seqw.size());
    out << std::endl;
    out << "Neff = " << hmm.neff;

    return out;
}

}  // namespace cs

#endif  // CS_PO_HMM_INL_H_
