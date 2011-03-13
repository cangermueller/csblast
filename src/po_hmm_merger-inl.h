// Copyright 2010, Andreas Biegert

#ifndef CS_PO_HMM_MERGER_INL_H_
#define CS_PO_HMM_MERGER_INL_H_

#include "po_hmm_merger.h"

namespace cs {

// Constructs a new merger for PO-HMMs 'x' and 'y'.
template<class Abc>
POHmmMerger<Abc>::POHmmMerger(const POHmm<Abc>& x, const POHmm<Abc>& y, double gapo, double gape)
        : x_(x),
          y_(y),
          MM_vertex_(x.size() + 1, y.size() + 1),
          MI_vertex_(x.size() + 1, y.size() + 1),
          IM_vertex_(x.size() + 1, y.size() + 1),
          gap_open(gapo),
          gap_extend(gape) {

    // New PO-HMM z should contain all sequences contained in x and y
    z_.seqs.insert(z_.seqs.end(), x_.seqs.begin(), x_.seqs.end());
    z_.seqs.insert(z_.seqs.end(), y_.seqs.begin(), y_.seqs.end());

    // Set gap open and gap extension probabilities of merged PO-HMM
    z_.gap_open   = gap_open;
    z_.gap_extend = gap_extend;

    // Add missing properties to START-END vertex
    for (size_t k = 0; k < z_.seqs.size(); ++k)
        z_.g[kStartEndVertex].col.push_back(std::make_pair(k, -1));

    // Add all vertices of PO-HMMs x to new PO-HMM z
    VertexIter vi, vertex_end;
    VertexIndex index;
    for (tie(vi, vertex_end) = vertices(x_.g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex)
            add_vertex(VertexProperties<Abc>(x_.g[*vi]), z_.g);
    }

    // Add all vertices of PO-HMMs y to new PO-HMM z
    for (tie(vi, vertex_end) = vertices(y_.g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex) {
            Vertex v = add_vertex(VertexProperties<Abc>(y_.g[*vi]), z_.g);
            for (ColIter ri = z_.g[v].col.begin(); ri != z_.g[v].col.end(); ++ri)
                ri->first += x_.seqs.size();
        }
    }

    // Add all edges of PO-HMMs x to new PO-HMM z
    EdgeIter ei, ei_end;
    for (tie(ei, ei_end) = edges(x_.g); ei != ei_end; ++ei) {
        size_t source_vertex = index[source(*ei, x_.g)];
        size_t target_vertex = index[target(*ei, x_.g)];
        add_edge(source_vertex, target_vertex, EdgeProperties(x_.g[*ei]), z_.g);
    }

    // Add all edges of PO-HMMs y to new PO-HMM z
    for (tie(ei, ei_end) = edges(y_.g); ei != ei_end; ++ei) {
        size_t source_vertex = index[source(*ei, y_.g)];
        size_t target_vertex = index[target(*ei, y_.g)];
        if (source_vertex != kStartEndVertex) source_vertex += x_.size();
        if (target_vertex != kStartEndVertex) target_vertex += x_.size();
        add_edge(source_vertex, target_vertex, EdgeProperties(y_.g[*ei]), z_.g);
    }
}

template<class Abc>
void POHmmMerger<Abc>::AddAlignment(const PairAlignment& ali) {
    PairAlignment full_ali(ali);                    // resulting ali with indels
    AlignmentPath& path = full_ali.path;            // quick access to path
    std::set<ColPair> res;                          // total residue set along path
    std::vector<int> last_r(x_.seqs.size(), -1);    // last residue for each seq
    std::vector<size_t> last_i(x_.seqs.size(), 0);  // and its vertex index
    VertexVec skipped;                              // skipped vertices in topo order
    std::vector<size_t> kvec;                       // keeps track of k indices
    std::vector<size_t>::iterator it;               // points to end of uniqe kvec interval

    LOG(ERROR) << "Adding alignment " << aliw_.size() << " to merger ...";
    LOG(INFO) << ali;
    LOG(ERROR) << "Extending alignment path for PO-HMM x ...";

    // Mark used vertices and determine initial residue set for PO-HMM x
    Bitset used_x(x_.size() + 1);
    used_x.set(kStartEndVertex);
    for (ConstPathIter s = ali.path.begin(); s != ali.path.end(); ++s) {
        if (s->pstate == MM || s->pstate == MI) {
            used_x.set(s->i);
            res.insert(x_.g[s->i].col.begin(), x_.g[s->i].col.end());
        }
    }

    // Construct an alignment path that fully covers PO-HMM x without skipping indels
    const Graph& x = x_.g;
    for (size_t s = 0; s < path.size(); ++s) {
        if (path[s].pstate == IM) continue;
        Vertex v = path[s].i;
        LOG(ERROR) << "s=" << s << " v=" << v << " path-size=" << path.size();

        // Check for missing residues in each sequence
        kvec.clear();
        for (ConstColIter r = x[v].col.begin(); r != x[v].col.end(); ++r) {
            if (last_r[r->first] + 1 != r->second) {
                kvec.push_back(last_i[r->first]);
                LOG(ERROR) << "Missing residues in sequence '" << x_.seqs[r->first].header() << "' of PO-HMM x:  k=" << kvec.back()
                           << "  v=" << x_.vertices[kvec.back()] << "  last_r=" << last_r[r->first] << "  curr_r=" << r->second;
            }
        }

        if (!kvec.empty()) {
            // One or more missing residues => recover them by adding skipped vertices
            std::sort(kvec.begin(), kvec.end());
            it = std::unique(kvec.begin(), kvec.end());
            kvec.resize(it - kvec.begin());
            skipped.clear();
            size_t v_end = v;
            for (std::vector<size_t>::reverse_iterator ki = kvec.rbegin(); ki != kvec.rend(); ++ki) {
                VertexVec skipped_k = GetSkippedVertices(x_.vertices[*ki], v_end, x_, used_x, res);
                for (typename VertexVec::iterator si = skipped_k.begin(); si != skipped_k.end(); ++si) used_x.set(*si);
                skipped.insert(skipped.begin(), skipped_k.begin(), skipped_k.end());
                v_end = x_.vertices[*ki];
            }

            // Find path index t of step to the right of vertex k and j index of step t-1
            size_t k = kvec.front();
            size_t t = s;
            size_t j = kStartEndVertex;
            while (t > 0) {
                if (x_.indices[path[t-1].i] == k && path[t-1].pstate != IM) {
                    j = path[t-1].j;
                    break;
                }
                t--;
            }
            LOG(INFO) << "t=" << t << " j=" << j;
            VertexVecIter vi = skipped.begin();
            while (!skipped.empty() && (vi != skipped.end() || t < s)) {
                if (vi != skipped.end() && x_.indices[*vi] < x_.indices[path[t].i]) {
                    LOG(INFO) << " inserting i=" << *vi << " j=" << j;
                    // Insert artifical MI step before step t
                    path.insert(path.begin() + t, Step(*vi, j, MI, 0.0));
                    used_x.set(*vi);
                    for (ConstColIter r = x[*vi].col.begin(); r != x[*vi].col.end(); ++r) {
                        LOG(INFO) << "  s=" << r->first << " r=" << r->second
                                  << " rlast=" << last_r[r->first] << " seq[s][r]="
                                  << Abc::kIntToChar[x_.seqs[r->first][r->second]];
                        last_r[r->first] = r->second;
                        last_i[r->first] = x_.indices[*vi];
                    }
                    vi++;
                    s++;
                    t++;
                } else {
                    LOG(INFO) << " skipping i=" << path[t].i << " j=" << path[t].j;
                    LOG(INFO) << "  t=" << t << " u=" << path[t].i;
                    Vertex u = path[t].i;
                    for (ConstColIter r = x[u].col.begin(); r != x[u].col.end(); ++r) {
                        LOG(INFO) << "  s=" << r->first << " r=" << r->second
                                  << " rlast=" << last_r[r->first] << " seq[s][r]="
                                  << Abc::kIntToChar[x_.seqs[r->first][r->second]];
                        last_r[r->first] = r->second;
                        last_i[r->first] = x_.indices[u];
                    }
                    j = path[t].j;
                    t++;  // insert all remaining skipped vertices after this one
                }
            }
        }

        // No missing residues => update 'last' vectors
        for (ConstColIter r = x[v].col.begin(); r != x[v].col.end(); ++r) {
            LOG(INFO) << " s=" << r->first << " r=" << r->second
                      << " rlast=" << last_r[r->first] << " seq[s][r]="
                      << Abc::kIntToChar[x_.seqs[r->first][r->second]];
            assert_eq(last_r[r->first] + 1, r->second);
            last_r[r->first] = r->second;
            last_i[r->first] = x_.indices[path[s].i];
        }
    }

    // Check for missing residues at the end of all sequences
    kvec.clear();
    for (size_t l = 0; l < x_.seqs.size(); ++l) {
        if (last_r[l] != static_cast<int>(x_.seqs[l].length() - 1)) {
            kvec.push_back(last_i[l]);
        }
    }

    if (!kvec.empty()) {
        // One or more one missing residues => recover them by adding skipped vertices
        std::sort(kvec.begin(), kvec.end());
        it = std::unique(kvec.begin(), kvec.end());
        kvec.resize(it - kvec.begin());
        skipped.clear();
        size_t v_end = kStartEndVertex;
        for (std::vector<size_t>::reverse_iterator ki = kvec.rbegin(); ki != kvec.rend(); ++ki) {
            VertexVec skipped_k = GetSkippedVertices(x_.vertices[*ki], v_end, x_, used_x, res);
            for (typename VertexVec::iterator si = skipped_k.begin(); si != skipped_k.end(); ++si) used_x.set(*si);
            skipped.insert(skipped.begin(), skipped_k.begin(), skipped_k.end());
            v_end = x_.vertices[*ki];
        }

        // Find path index t of step to the right of vertex k and j index of step t-1
        size_t k = kvec.front();
        size_t t = path.size();
        size_t j = kStartEndVertex;
        while (t > 0) {
            if (x_.indices[path[t-1].i] == k && path[t-1].pstate != IM) {
                j = path[t-1].j;
                break;
            }
            t--;
        }
        LOG(INFO) << "t=" << t << " j=" << j;

        VertexVecIter vi = skipped.begin();
        while (!skipped.empty() && (vi != skipped.end() || t < path.size())) {
            LOG(INFO) << "t=" << t;
            if (vi != skipped.end() &&
                (t == path.size() || x_.indices[*vi] < x_.indices[path[t].i])) {
                LOG(INFO) << " inserting i=" << *vi << " j=" << j;
                // Insert artifical MI step before step t
                path.insert(path.begin() + t, Step(*vi, j, MI, 0.0));
                used_x.set(*vi);
                for (ConstColIter r = x[*vi].col.begin(); r != x[*vi].col.end(); ++r) {
                    last_r[r->first] = r->second;
                    last_i[r->first] = x_.indices[*vi];
                }
                vi++;
                t++;
            } else {
                LOG(INFO) << " skipping i=" << path[t].i << " j=" << path[t].j;
                Vertex u = path[t].i;
                for (ConstColIter r = x[u].col.begin(); r != x[u].col.end(); ++r) {
                    last_r[r->first] = r->second;
                    last_i[r->first] = x_.indices[u];
                }
                j = path[t].j;
                t++;  // insert all remaining skipped vertices after this one
            }
        }
    }

    LOG(INFO) << full_ali;
    LOG(ERROR) << "Extending alignment path for PO-HMM y ...";

    // Reuse variables for extending alignment path with respect to PO-HMM y
    res.clear();
    last_r = std::vector<int>(y_.seqs.size(), -1);
    last_i = std::vector<size_t>(y_.seqs.size(), 0);

    // Mark used vertices and determine initial residue set for PO-HMM y
    Bitset used_y(y_.size() + 1);
    used_y.set(kStartEndVertex);
    for (ConstPathIter s = ali.path.begin(); s != ali.path.end(); ++s) {
        if (s->pstate == MM || s->pstate == IM) {
            used_y.set(s->j);
            res.insert(y_.g[s->j].col.begin(), y_.g[s->j].col.end());
        }
    }

    // Construct an alignment path that fully covers PO-HMM y without skipping indels
    const Graph& y = y_.g;
    for (size_t s = 0; s < path.size(); ++s) {
        if (path[s].pstate == MI) continue;
        Vertex v = path[s].j;
        LOG(ERROR) << "s=" << s << " v=" << v << " path-size=" << path.size();

        // Check for missing residues in each sequence
        kvec.clear();
        for (ConstColIter r = y[v].col.begin(); r != y[v].col.end(); ++r) {
            if (last_r[r->first] + 1 != r->second) {
                kvec.push_back(last_i[r->first]);
                LOG(ERROR) << "Missing residues in sequence '" << y_.seqs[r->first].header() << "'  of PO-HMM y:  k="  << kvec.back()
                           << "  v=" << y_.vertices[kvec.back()] << "  last_r=" << last_r[r->first] << "  curr_r=" << r->second;
            }
        }

        if (!kvec.empty()) {
            // One or more missing residues => recover them by adding skipped vertices
            std::sort(kvec.begin(), kvec.end());
            it = std::unique(kvec.begin(), kvec.end());
            kvec.resize(it - kvec.begin());
            skipped.clear();
            size_t v_end = v;
            for (std::vector<size_t>::reverse_iterator ki = kvec.rbegin(); ki != kvec.rend(); ++ki) {
                VertexVec skipped_k = GetSkippedVertices(y_.vertices[*ki], v_end, y_, used_y, res);
                for (typename VertexVec::iterator si = skipped_k.begin(); si != skipped_k.end(); ++si) used_y.set(*si);
                skipped.insert(skipped.begin(), skipped_k.begin(), skipped_k.end());
                v_end = y_.vertices[*ki];
            }

            // Find path index t of step to the right of vertex k and i index of step t-1
            size_t k = kvec.front();
            size_t t = s;
            size_t i = kStartEndVertex;
            while (t > 0) {
                if (y_.indices[path[t-1].j] == k && path[t-1].pstate != MI) {
                    i = path[t-1].i;
                    break;
                }
                t--;
            }
            LOG(INFO) << "t=" << t << " i=" << i;
            VertexVecIter vj = skipped.begin();
            while (!skipped.empty() && (vj != skipped.end() || t < s)) {
                if (vj != skipped.end() && y_.indices[*vj] < y_.indices[path[t].j]) {
                    LOG(INFO) << " inserting i=" << i << " j=" << *vj;
                    // Insert artifical IM step before step t
                    path.insert(path.begin() + t, Step(i, *vj, IM, 0.0));
                    used_y.set(*vj);
                    for (ConstColIter r = y[*vj].col.begin(); r != y[*vj].col.end(); ++r) {
                        last_r[r->first] = r->second;
                        last_i[r->first] = y_.indices[*vj];
                    }
                    vj++;
                    s++;
                    t++;
                } else {
                    LOG(INFO) << " skipping i=" << path[t].i << " j=" << path[t].j;
                    Vertex u = path[t].j;
                    for (ConstColIter r = y[u].col.begin(); r != y[u].col.end(); ++r) {
                        last_r[r->first] = r->second;
                        last_i[r->first] = y_.indices[u];
                    }
                    i = path[t].i;
                    t++;  // insert all remaining skipped vertices after this one
                }
            }
        }

        // No missing residues => update 'last' vectors
        for (ConstColIter r = y[v].col.begin(); r != y[v].col.end(); ++r) {
            LOG(INFO) << " s=" << r->first << " r=" << r->second
                      << " rlast=" << last_r[r->first];
            assert_eq(last_r[r->first] + 1, r->second);
            last_r[r->first] = r->second;
            last_i[r->first] = y_.indices[path[s].j];
        }
    }

    // Check for missing residues at the end of all sequences
    kvec.clear();
    for (size_t l = 0; l < y_.seqs.size(); ++l) {
        if (last_r[l] != static_cast<int>(y_.seqs[l].length() - 1)) {
            kvec.push_back(last_i[l]);
        }
    }

    if (!kvec.empty()) {
        // One or more one missing residues => recover them by adding skipped vertices
        std::sort(kvec.begin(), kvec.end());
        it = std::unique(kvec.begin(), kvec.end());
        kvec.resize(it - kvec.begin());
        skipped.clear();
        size_t v_end = kStartEndVertex;
        for (std::vector<size_t>::reverse_iterator ki = kvec.rbegin(); ki != kvec.rend(); ++ki) {
            VertexVec skipped_k = GetSkippedVertices(y_.vertices[*ki], v_end, y_, used_y, res);
            for (typename VertexVec::iterator si = skipped_k.begin(); si != skipped_k.end(); ++si) used_y.set(*si);
            skipped.insert(skipped.begin(), skipped_k.begin(), skipped_k.end());
            v_end = y_.vertices[*ki];
        }

        // Find path index t of step to the right of vertex k and j index of step t-1
        size_t k = kvec.front();
        size_t t = path.size();
        size_t i = kStartEndVertex;
        while (t > 0) {
            if (y_.indices[path[t-1].j] == k && path[t-1].pstate != MI) {
                i = path[t-1].i;
                break;
            }
            t--;
        }
        LOG(INFO) << "t=" << t << " i=" << i;

        VertexVecIter vj = skipped.begin();
        while (!skipped.empty() && (vj != skipped.end() || t < path.size())) {
            LOG(INFO) << "t=" << t;
            if (vj != skipped.end() &&
                (t == path.size() || y_.indices[*vj] < y_.indices[path[t].j])) {
                LOG(INFO) << " inserting i=" << i << " j=" << *vj;
                // Insert artifical IM step before step t
                path.insert(path.begin() + t, Step(i, *vj, IM, 0.0));
                used_y.set(*vj);
                for (ConstColIter r = y[*vj].col.begin(); r != y[*vj].col.end(); ++r) {
                    last_r[r->first] = r->second;
                    last_i[r->first] = y_.indices[*vj];
                }
                vj++;
                t++;
            } else {
                LOG(INFO) << " skipping i=" << path[t].i << " j=" << path[t].j;
                Vertex u = path[t].j;
                for (ConstColIter r = y[u].col.begin(); r != y[u].col.end(); ++r) {
                    last_r[r->first] = r->second;
                    last_i[r->first] = y_.indices[u];
                }
                i = path[t].i;
                t++;  // insert all remaining skipped vertices after this one
            }
        }
    }

    LOG(INFO) << full_ali;

    // Use the full alignment path as guide for merging vertices in x and y
    z_.g[kStartEndVertex].alis.set(aliw_.size());  // mark ali in START-END vertex
    for (ConstPathIter s = full_ali.path.begin(); s != full_ali.path.end(); ++s) {
        AddVertex(s->i, s->j, s->pstate);
    }
    z_.num_alis++;
    aliw_.push_back(ali.pmin);
}


template<class Abc>
typename POHmmMerger<Abc>::VertexVec POHmmMerger<Abc>::GetSkippedVertices(Vertex v1,
                                                                          Vertex v2,
                                                                          const POHmm<Abc>& hmm,
                                                                          const Bitset& path,
                                                                          const std::set<ColPair>& res) const {
    // Constants used throughout indel detection algorithm
    const Graph& g = hmm.g;
    const VertexVec& vertices = hmm.vertices;
    const IndexVec& indices = hmm.indices;
    size_t beg = indices[v1];
    size_t end = (v2 == kStartEndVertex) ? vertices.size() : indices[v2];
    VertexVec skipped_vertices; // return value

    LOG(ERROR) << "Finding skipped vertices between v1=" << v1 << " and v2=" << v2;

    // Construct bitmask by AND gating alis of used vertices between 'beg' and 'end'
    AliBitset mask = g[v1].alis & g[v2].alis;
    LOG(ERROR) << "mask=" << mask.to_string();
    for (size_t k = beg + 1; k < end; ++k) {
        if (path.test(vertices[k])) {
            mask &= g[vertices[k]].alis;
            LOG(ERROR) << "v=" << vertices[k] << " mask = " << mask.to_string();
        }
    }

    // Choose mask as AND gating of 'beg' and 'end' alis in case no continuous path is found
    if (mask.none()) mask = g[v1].alis & g[v2].alis;
    LOG(ERROR) << "mask=" << mask.to_string();

    // Forward pass from left to right
    Bitset off_fwd(num_vertices(g));
    std::set<Vertex> set_fwd;
    for (size_t k = beg + 1; k < end; ++k) {
        Vertex v = vertices[k];
        if ((g[v].alis & mask).none()) { off_fwd.set(v); continue; }
        // Check residue clashes
        bool found = false;
        for (ConstColIter r = g[v].col.begin(); r != g[v].col.end(); ++r)
            if (res.find(*r) != res.end()) { found = true; break; }
        if (found) { off_fwd.set(v); continue; }
        // Check edge connectivity
        bool accepted = false;
        InEdgeIter ei, edge_end;
        for (tie(ei, edge_end) = in_edges(v, g); ei != edge_end; ++ei) {
            Vertex u = source(*ei, g);
            if (path.test(u) || (indices[u] > beg && !off_fwd.test(u))) {
                accepted = true;
                break;
            }
        }
        if (accepted) set_fwd.insert(v);
        else off_fwd.set(v);
    }

    // Backward pass from right to left
    Bitset off_bwd(num_vertices(g));
    std::set<Vertex> set_bwd;
    for (size_t k = end - 1; k > beg; --k) {
        Vertex v = vertices[k];
        if ((g[v].alis & mask).none()) { off_bwd.set(v); continue; }
        // Check residue clashes
        bool found = false;
        for (ConstColIter r = g[v].col.begin(); r != g[v].col.end(); ++r)
            if (res.find(*r) != res.end()) { found = true; break; }
        if (found) { off_bwd.set(v); continue; }
        // Check edge connectivity
        bool accepted = false;
        OutEdgeIter ei, edge_end;
        for (tie(ei, edge_end) = out_edges(v, g); ei != edge_end; ++ei) {
            Vertex u = target(*ei, g);
            if (end != kStartEndVertex) {
                if (path.test(u) || (indices[u] < end && !off_bwd.test(u))) {
                    accepted = true;
                    break;
                }
            } else {
                if (path.test(u) || !off_bwd.test(u)) {
                    accepted = true;
                    break;
                }
            }
        }
        if (accepted) set_bwd.insert(v);
        else off_bwd.set(v);
    }
    LOG(ERROR) << "set_fwd = " << StringifyContainer(set_fwd);
    LOG(ERROR) << "set_bwd = " << StringifyContainer(set_bwd);

    // Return empty vector in case no matching vertices were found
    if (set_fwd.empty() || set_bwd.empty()) return skipped_vertices;

    // Determine set of candidate vertices as intersection of fwd and bwd sets
    std::set<Vertex> cand;
    set_intersection(set_fwd.begin(), set_fwd.end(), set_bwd.begin(), set_bwd.end(), inserter(cand, cand.begin()));
    LOG(ERROR) << "cand = " << StringifyContainer(cand);

    // Return empty vector in case there are no candidates
    if (cand.empty()) return skipped_vertices;

    // Determine position of lowest set bit in 'mask'
    size_t lowest_pos = 0;
    for (size_t pos = 0; pos < kMaxAliCapacity; ++pos) {
        if (mask.test(pos)) {
            lowest_pos = pos;
            break;
        }
    }
    LOG(ERROR) << "lowest_pos = " << lowest_pos;

    // // Build tuples vector of candidate vertices
    // std::vector<CountAlisVertexTuple> tuples;
    // for (size_t k = beg + 1; k < end; ++k) {
    //     Vertex v = vertices[k];
    //     if (cand.find(v) != cand.end())
    //         tuples.push_back(boost::make_tuple(g[v].alis.count(), g[v].alis.to_string(), v));
    // }
    // std::sort(tuples.begin(), tuples.end());  // sort vertices by alignment bitcount

    // // Find characteristic bitset of path with maximum flow
    // AliBitset alis_max, alis_all;
    // double flow_max = 0.0;
    // for (size_t i = 0; i < tuples.size() && alis_all != mask; ++i) {
    //     Vertex v = boost::get<2>(tuples[i]);
    //     if ((g[v].alis & alis_max).none() && (g[v].alis & alis_all).none() && g[v].flow > flow_max) {
    //         flow_max = g[v].flow;
    //         alis_max = g[v].alis;
    //         LOG(ERROR) << strprintf("v=%-5zu  flow_max=%6.4f  alis_max=%s", v, flow_max, alis_max.to_string().c_str());
    //     }
    //     alis_all |= g[v].alis;
    // }
    // LOG(ERROR) << "alis_max=" << alis_max.to_string();

    // // Extract vertices of path with maximum flow
    // for (size_t k = beg + 1; k < end; ++k) {
    //     Vertex v = vertices[k];
    //  if (cand.find(v) != cand.end() && (g[v].alis & alis_max).any() && !path.test(v))
    //      // if (cand.find(v) != cand.end() && (g[v].alis & alis_max) == alis_max && !path.test(v))
    //         skipped_vertices.push_back(v);
    // }

    // Extract vertices with set 'lowest_pos' alignment bit
    for (size_t k = beg + 1; k < end; ++k) {
        Vertex v = vertices[k];
        if (cand.find(v) != cand.end() && g[v].alis.test(lowest_pos) && !path.test(v))
            skipped_vertices.push_back(v);
    }

    LOG(ERROR) << "Skipped vertices:";
    LOG(ERROR) << StringifyContainer(skipped_vertices);
    return skipped_vertices;
}

template<class Abc>
void POHmmMerger<Abc>::AddVertex(size_t i, size_t j, uint8_t pstate) {
    InEdgeIter ie, ie_end;
    OutEdgeIter oe, oe_end;
    VertexIndex index;

    // Create vertex anew or fetch it from cache
    Vertex v = FindOrCreateVertex(i, j, pstate);

    // Set alignment-bit to keep track which alignments used this vertex
    z_.g[v].alis.set(aliw_.size());

    // Build bitset of v's in-edge and out-edge  connections for quick checking
    // TODO: could we make this a member instead of recomputing it every time?
    Bitset v_in(z_.size() + 1);
    for (tie(ie, ie_end) = in_edges(v, z_.g); ie != ie_end; ++ie)
        v_in.set(index[source(*ie, z_.g)]);
    Bitset v_out(z_.size() + 1);
    for (tie(oe, oe_end) = out_edges(v, z_.g); oe != oe_end; ++oe)
        v_out.set(index[target(*oe, z_.g)]);

    // Add in-edges and out-edges of vertex i to vertex v
    if (pstate == MM || pstate == MI) {
        // Copy in-edges of vertex i to vertex v
        for (tie(ie, ie_end) = in_edges(i, z_.g); ie != ie_end; ++ie) {
            size_t k = index[source(*ie, z_.g)];
            if (!v_in.test(k)) {
                add_edge(k, v, EdgeProperties(), z_.g);
                v_in.set(k);
            }
        }
        // Copy out-edges of vertex i to vertex v
        for (tie(oe, oe_end) = out_edges(i, z_.g); oe != oe_end; ++oe) {
            size_t l = index[target(*oe, z_.g)];
            if (!v_out.test(l)) {
                add_edge(v, l, EdgeProperties(), z_.g);
                v_out.set(l);
            }
        }
    }

    // Add in-edges and out-edges of vertex j to vertex v
    if (pstate == MM || pstate == IM) {
        // Copy in-edges of vertex j to vertex v
        for (tie(ie, ie_end) = in_edges(j + x_.size(), z_.g); ie != ie_end; ++ie) {
            size_t k = index[source(*ie, z_.g)];
            if (!v_in.test(k)) {
                add_edge(k, v, EdgeProperties(), z_.g);
                v_in.set(k);
            }
        }
        // Copy out-edges of vertex j to vertex v
        for (tie(oe, oe_end) = out_edges(j + x_.size(), z_.g); oe != oe_end; ++oe) {
            size_t l = index[target(*oe, z_.g)];
            if (!v_out.test(l)) {
                add_edge(v, l, EdgeProperties(), z_.g);
                v_out.set(l);
            }
        }
    }
}

template<class Abc>
size_t POHmmMerger<Abc>::FindOrCreateVertex(size_t i, size_t j, uint8_t pstate) {
    const size_t x_nseqs = x_.seqs.size();
    const VertexProperties<Abc>& vi = x_.g[i];
    const VertexProperties<Abc>& vj = y_.g[j];

    // Check if vertex has already been constructed
    switch (pstate) {
        case MM:
            if (MM_vertex_.test(i,j)) return MM_vertex_.get(i,j);
            break;
        case MI:
            if (MI_vertex_.test(i,j)) return MI_vertex_.get(i,j);
            break;
        case IM:
            if (IM_vertex_.test(i,j)) return IM_vertex_.get(i,j);
            break;
    }

    // Vertex not yet contained in z, we contruct it by merging vertices i and j
    Vertex v = add_vertex(VertexProperties<Abc>(), z_.g);

    // Set pair state and aligned positions so we can assign proper phylogeny based probs later on
    z_.g[v].pstate = pstate;
    z_.g[v].i = i;
    z_.g[v].j = j;

    // Merge col entries
    if (pstate == MM || pstate == MI) {
        for (ConstColIter ri = vi.col.begin(); ri != vi.col.end(); ++ri)
            z_.g[v].col.push_back(std::make_pair(ri->first, ri->second));
        z_.g[v].skip_count = pstate == MM ? 0 : vi.skip_count + 1;
    }
    if (pstate == MM || pstate == IM) {
        for (ConstColIter ri = vj.col.begin(); ri != vj.col.end(); ++ri)
            z_.g[v].col.push_back(std::make_pair(ri->first + x_nseqs, ri->second));
        z_.g[v].skip_count = pstate == MM ? 0 : vj.skip_count + 1;
    }

    // Store new vertex so we can reuse it later on
    switch (pstate) {
        case MM:
            MM_vertex_.set(i, j, v);
            break;
        case MI:
            MI_vertex_.set(i, j, v);
            break;
        case IM:
            IM_vertex_.set(i, j, v);
            break;
    }
    return v;
}

template<class Abc>
POHmm<Abc> POHmmMerger<Abc>::Finalize() {
    const size_t xy_offset = x_.size() + y_.size();
    POHmm<Abc> z(z_);  // z is the finalized copy of 'z_' without dummy vertices
    VertexIter vi, vertex_end;
    EdgeIter ei, ei_end;
    VertexIndex index;

    // Clear z and add START-END vertex before copying "upper" vertices of 'z_'
    z.g = Graph();
    add_vertex(VertexProperties<Abc>(), z.g);
    z.g[kStartEndVertex] = z_.g[kStartEndVertex];

    // Add "upper" vertices of PO-HMMs 'z_' to 'z' and set flow
    for (tie(vi, vertex_end) = vertices(z_.g); vi != vertex_end; ++vi) {
        if (index[*vi] != kStartEndVertex && index[*vi] > xy_offset)
            add_vertex(VertexProperties<Abc>(z_.g[*vi]), z.g);
    }

    // Add all edges between upper vertices in 'z_' to 'z'
    for (tie(ei, ei_end) = edges(z_.g); ei != ei_end; ++ei) {
        size_t u = index[source(*ei, z_.g)];
        size_t v = index[target(*ei, z_.g)];
        if ((u > xy_offset || u == kStartEndVertex) &&
            (v > xy_offset || v == kStartEndVertex) &&
            (z_.g[u].alis & z_.g[v].alis).any()) {
            size_t uu = (u == kStartEndVertex) ? kStartEndVertex : u - xy_offset;
            size_t vv = (v == kStartEndVertex) ? kStartEndVertex : v - xy_offset;
            add_edge(uu, vv, EdgeProperties(z_.g[*ei]), z.g);
        }
    }

    // Normalize alignment weights and assign to z
    z.aliw = aliw_;
    Normalize(&z.aliw[0], z.aliw.size());
    z.num_alis = aliw_.size();

    // Calculate vertex flows by summing weights of alignments at each vertex
    for (tie(vi, vertex_end) = vertices(z.g); vi != vertex_end; ++vi) {
        z.g[*vi].flow = 0.0;
        for (size_t n = 0; n < z.num_alis; ++n)
            if (z.g[*vi].alis.test(n)) z.g[*vi].flow += z.aliw[n];
    }

    // Set alignment bits at START-END vertex
    z.g[kStartEndVertex].alis.reset();
    for (size_t a = 0; a < z.num_alis; ++a) z.g[kStartEndVertex].alis.set(a);

    // Use vertex flows to calculate global seq-weights needed for transitions probs
    z.Init();

    return z;
}

}  // namespace cs

#endif  // CS_PO_HMM_MERGER_INL_H_
