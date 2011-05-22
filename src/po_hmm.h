// Copyright 2009, Andreas Biegert

#ifndef CS_PO_HMM_H_
#define CS_PO_HMM_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>

#include "alignment-inl.h"
#include "profile_column.h"
#include "sequence-inl.h"
#include "sparse_profile.h"

using boost::adjacency_list;
using boost::vecS;
using boost::bidirectionalS;
using boost::graph_traits;
using boost::filtered_graph;
using boost::property_map;
using boost::vertex_index_t;
using boost::vertex_index;

namespace cs {

// Forward declarations
template<class Abc>
class POHmm;

template<class Abc>
class ContextLibrary;

template<class Abc>
class Emission;

// Index of start vertex
static const size_t kStartEndVertex = 0;
// Maximial number of alignment paths that can be stored in a PO-HMM
static const size_t kMaxAliCapacity = 64;
// Maximial number of sequences that can be stored in a PO-HMM
static const size_t kMaxSeqCapacity = 64;
// default gap-open probability
static const double kGapOpen = 0.025;
// Default gap-extension probability
static const double kGapExt  = 0.75;

// Typedefs needed throughout all code dealing with PO-HMMs
typedef std::bitset<kMaxAliCapacity> AliBitset;  // keeps track of alignment paths
typedef std::pair<int, int> ColPair;             // sequence-residue pair
typedef std::vector<ColPair> ColVec;             // aligned residues at vertex

// Iterators
typedef ColVec::iterator ColIter;
typedef ColVec::const_iterator ConstColIter;

// Pair states needed for backtracing
enum PairState { MM=0, MI=1, IM=2 };

// Strategy class for initializing a PO-HMM
template<class Abc>
class POHmmInit {
  public:
    POHmmInit() {}
    virtual ~POHmmInit() {}
    virtual void operator() (POHmm<Abc>& hmm) const = 0;
};

// POD class representing a POG vertex which holds among others information about
// transition probs (M2M, M2I, etc.), aligned residues, and abstract state probs.
template<class Abc>
struct VertexProperties {
    // Constructs an empty vertex
    VertexProperties()
            : M2M(0.0),
              M2I(0.0),
              I2M(0.0),
              I2I(0.0),
              counts(0.0),
              probs(0.0),
              neff(0),
              flow(0.0),
              i(0),
              j(0),
              pstate(0),
              skip_count(0) {}

    double M2M, M2I, I2M, I2I;     // transition probs
    SparseProfileCol contexts;     // sparse list of context probabilities
    ProfileColumn<Abc> counts;     // unnormalized residue counts based on seq weights WITHOUT pseudocounts
    ProfileColumn<Abc> probs;      // same as above but normalized to one and possibly WITH pseudocounts
    ColVec col;                    // aligned sequences residues in this column
    double neff;                   // number of effective sequences
    double flow;                   // weights sum of alignments that use this vertex
    AliBitset alis;                // for flagging alis that use this vertex
    size_t i, j;                   // aligned positions in merged HMM
    uint8_t pstate;                // pair state of merged vertex
    size_t skip_count;             // counts how often this vertex has been skipped over in progressive alignment
};

// POD class representing a POG edge
struct EdgeProperties {
    EdgeProperties(double w = 1.0, double wr = 1.0) : weight(w), weight_rev(wr) {}
    double weight,  weight_rev;
};


// A container class representing a partial-order HMM.
template<class Abc>
class POHmm {
  public:
    // Typedefs
    typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProperties<Abc>,
                           EdgeProperties> Graph; // this is our graph!
    typedef graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor Vertex;
    typedef typename GraphTraits::edge_descriptor Edge;
    typedef typename GraphTraits::vertex_iterator VertexIter;
    typedef typename GraphTraits::edge_iterator EdgeIter;
    typedef typename GraphTraits::in_edge_iterator InEdgeIter;
    typedef typename GraphTraits::out_edge_iterator OutEdgeIter;
    typedef typename property_map<Graph, vertex_index_t>::type VertexIndex;
    typedef typename std::vector<Sequence<Abc> > SeqVec;
    typedef typename std::vector<Vertex> VertexVec;
    typedef typename VertexVec::iterator VertexVecIter;
    typedef typename VertexVec::const_iterator ConstVertexVecIter;
    typedef typename std::vector<size_t> IndexVec;
    typedef typename std::vector<double> WeightsVec;

    // Constructs an empty PO-HMM containing only the START-END vertex.
    POHmm(double gapo = kGapOpen, double gape = kGapExt);

    // Constructs a PO-HMM from single sequence
    POHmm(const Sequence<Abc>& seq, double gapo = kGapOpen, double gape = kGapExt);

    // Constructs a PO-HMM from a multiple alignment consisting of match columns only.
    POHmm(const Alignment<Abc>& ali, double gapo = kGapOpen, double gape = kGapExt);

    // Returns number of vertices in the POG not counting the begin-end vertex
    size_t size() const { return num_vertices(g) - 1; }

    // Calculates profiles of relative and absolute counts based on sequence
    void AssignCountsAndProbs();

    // Assigns transition probs based on 'gap_open', 'gap_extend', and gap patterns
    void AssignTransitionProbs();

    // Calculates edge probs based on branching probs and global sequence weights.
    void AssignEdgeWeights();

    // For each vertex calculates posterior probs of context states in 'lib' and
    // stores all probs above 'cutoff' in sparse profile column.
    void AssignContextStateProbs(const ContextLibrary<Abc>& lib, const Emission<Abc>& emission, double cutoff = 0.01);

    // Calculates and assigns probs of ancestral residues for each vertex.
    void AssignAncestralProbs(const POHmm<Abc>& x,
                              const POHmm<Abc>& y,
                              const SubstitutionMatrix<Abc>& matx,
                              const SubstitutionMatrix<Abc>& maty);

    // Sets emission probs of vertices with skip-count greater than 'max_skip_count' to zero.
    void DisableSkipped(size_t max_skip_count);

    // Calculates global seq-weights using the flow-weighted maximum-entropy method.
    void CalculateSequenceWeights();

    // Returns the n-th suboptimal PO-HMM alignment as multiple alignment of all seqs
    Alignment<Abc> GetAlignment(size_t n = 0) const;

    // Initializes seq-weights, edge weights, transitions, counts, and probs
    // based on 'col' entries, alignment weights, and flows
    void Init();

    // Returns an iterator pointing to beginning of topologically sorted vertices.
    VertexVecIter begin() { return vertices.begin(); }

    // Returns an iterator pointing past the end of topologically sorted vertices.
    VertexVecIter end() { return vertices.end(); }

    // Returns a const iter pointing to beginning of topologically sorted vertices.
    ConstVertexVecIter begin() const { return vertices.begin(); }

    // Returns a const iter pointing past the end of topologically sorted vertices.
    ConstVertexVecIter end() const { return vertices.end(); }

    // Checks if there are no additional or missing residues in PO-HMM alignments
    // with respect to original input sequences.
    void SanityCheck() const;

    Graph g;             // this is the actual graph data structure
    SeqVec seqs;         // sequences already contained in POG
    WeightsVec seqw;     // global sequence weights
    WeightsVec aliw;     // normalized pmin's of all included alignments
    size_t num_alis;     // number of pairwise alignments encoded in POG
    double neff;         // overall Neff
    double gap_open;     // default gap-open probability
    double gap_extend;   // default gap-extension probability
    VertexVec vertices;  // topologically sorted vertices
    IndexVec indices;    // maps each vertex to its index in topological order

  private:
    // Functor for creating a filtered graph without transitions into end-state.
    // This is needed because topological sort requires a DAG without cycles!
    struct EdgeFilter {
        EdgeFilter() : g(NULL) { }
        EdgeFilter(const Graph* graph_ptr) : g(graph_ptr) { }

        bool operator()(const Edge& e) const {
            VertexIndex index;
            return index[target(e, *g)] != kStartEndVertex;
        }

        const Graph* g;
    };

    struct InEdgeIterCompare :
            public std::binary_function<InEdgeIter, InEdgeIter, bool> {
        InEdgeIterCompare(const Graph& graph, const Vector<int>& vertex2index)
                : g(graph), v2idx(vertex2index) {}

        bool operator()  (const InEdgeIter& ei, const InEdgeIter& ej) const {
            return (v2idx[source(*ei, g)] > v2idx[source(*ej, g)]);
        }

        const Graph& g;
        const Vector<int>& v2idx;
    };

    struct OutEdgeIterCompare :
            public std::binary_function<OutEdgeIter, OutEdgeIter, bool> {
        OutEdgeIterCompare(const Graph& graph, const Vector<int>& vertex2index)
                : g(graph), v2idx(vertex2index) {}

        bool operator()  (const OutEdgeIter& ei, const OutEdgeIter& ej) const {
            return (v2idx[target(*ei, g)] < v2idx[target(*ej, g)]);
        }

        const Graph& g;
        const Vector<int>& v2idx;
    };
};


// Prints a PO-HMM in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const POHmm<Abc>& hmm);


}  // namespace cs

#endif  // CS_PO_HMM_H_
