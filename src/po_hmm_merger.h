// Copyright 2010, Andreas Biegert

#ifndef CS_PO_HMM_MERGER_H_
#define CS_PO_HMM_MERGER_H_

#include "po_hmm-inl.h"
#include "pairwise_aligner-inl.h"

namespace cs {

// Strategy class for merging two PO-HMMs into a new one using one or more
// pairwise alignments as guides.
template<class Abc>
class POHmmMerger {
  public:
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::Vertex Vertex;
    typedef typename POHmm<Abc>::Edge Edge;
    typedef typename POHmm<Abc>::VertexIter VertexIter;
    typedef typename POHmm<Abc>::VertexIndex VertexIndex;
    typedef typename POHmm<Abc>::EdgeIter EdgeIter;
    typedef typename POHmm<Abc>::OutEdgeIter OutEdgeIter;
    typedef typename POHmm<Abc>::InEdgeIter InEdgeIter;
    typedef typename POHmm<Abc>::VertexVec VertexVec;
    typedef typename POHmm<Abc>::IndexVec IndexVec;
    typedef typename POHmm<Abc>::VertexVecIter VertexVecIter;
    typedef boost::tuple<size_t, std::string, Vertex> CountAlisVertexTuple;

    // Constructs a new merger for PO-HMMs 'x' and 'y'.
    POHmmMerger(const POHmm<Abc>& x, const POHmm<Abc>& y, double gapo = kGapOpen, double gape = kGapExt);

    // Merges vertices according to given pairwise alignments.
    void AddAlignment(const PairAlignment& ali);

    // Indicates merger that no further alignments will be added and triggers
    // postprocessing of final PO-HMM. The resulting PO-HMM is returned.
    POHmm<Abc> Finalize();

    // Returns size of merged PO-HMM if we would finalize at current stage.
    size_t size() const { return z_.size() - x_.size() - y_.size(); }

  private:
    // Fetches given vertex from cache or constructs it by merging vertices i and j
    size_t FindOrCreateVertex(size_t i, size_t j, uint8_t pstate);

    // Incorporates vertex corresponding to given alignment step into PO-HMM z
    void AddVertex(size_t i, size_t j, uint8_t pstate);

    // Returns vertices that lie topologically between vertices i and j in graph g
    VertexVec GetSkippedVertices(Vertex v1,
                                 Vertex v2,
                                 const POHmm<Abc>& hmm,
                                 const Bitset& path,
                                 const std::set<ColPair>& res) const;

    const POHmm<Abc>& x_, y_;         // PO-HMMs that should be merged
    POHmm<Abc> z_;                    // merged PO-HMM as it is constructed
    SparseMatrix<Vertex> MM_vertex_;  // keep track of constructed MM vertices
    SparseMatrix<Vertex> MI_vertex_;  // keep track of constructed MI vertices
    SparseMatrix<Vertex> IM_vertex_;  // keep track of constructed IM vertices
    std::vector<double> aliw_;        // pmin's of included alignments
    double gap_open;                  // gap-open probability
    double gap_extend;                // gap-extension probability
};

}  // namespace cs

#endif  // CS_PO_HMM_MERGER_H_
