// Copyright 2010, Andreas Biegert

#ifndef CS_PAIRWISE_ALIGNER_H_
#define CS_PAIRWISE_ALIGNER_H_

#include "context_library-inl.h"
#include "po_hmm-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Many of these make up an alignment path.
struct Step {
    Step(size_t k, size_t l, uint8_t s, double p) : i(k), j(l), pstate(s), prob(p) {}
    size_t i, j;     // aligned vertices in PO-HMMs x and y
    uint8_t pstate;  // pair-state type
    double prob;     // posterior probability
};

// Keeps track of path through dynamic programing matrix
struct BacktraceCell {
    BacktraceCell(size_t i = 0, size_t j = 0, uint8_t ps = MM)
            : from_i(i), from_j(j), pstate(ps) {}
    size_t from_i, from_j;  // source vertices i and j in PO-HMMs x and y
    uint8_t pstate;         // pair-state type
};

// Keeps track of forward/backward probs
struct ProbCell {
    ProbCell() : MM(0.0), MI(0.0), IM(0.0) {}
    double MM, MI, IM;  // probs for being in MM, MI, and IM state
};

// Keeps track of score sum, average score, and path length during MAC alignment
struct ScoreCell {
    ScoreCell() : sum(0.0), avg(0.0), len(0) {}
    double sum, avg;
    size_t len;
};

// Typedefs used throughout all alignment algorithms
typedef Matrix<BacktraceCell> BacktraceMatrix;
typedef Matrix<ProbCell> ProbMatrix;
typedef Matrix<ScoreCell> ScoreMatrix;
typedef std::vector<Step> AlignmentPath;
typedef AlignmentPath::const_iterator ConstPathIter;
typedef AlignmentPath::iterator PathIter;
typedef std::pair<size_t, uint8_t> IndexStatePair;
typedef std::vector<IndexStatePair> IndexPath;

// Simple struct that holds all data from a pairwise alignment
struct PairAlignment {
    PairAlignment()
            : subopt(0), num_matches(0), sum_probs(0.0), avg_probs(0.0), pmin(1.0) {}
    size_t subopt;       // is incremented for suboptimal alignments, zero-based
    size_t num_matches;  // number of matches in the alignment
    double sum_probs;    // sum of posterior probabilities along path (excl. gapf!)
    double avg_probs;    // average posterior probabilities along path
    double pmin;         // minimum of posterior probabilities along path
    AlignmentPath path;  // path of aligned vertices
};

// POD that holds the two PO-HMMs to be aligned together with the dynamic
// programming matrices needed for forward, backward, and MAC alignment.
// Note that we copy both PO-HMMs since the aligner adds pseudocounts to their
// profiles and includes the null model context probabilities in the second PO-HMM.
template<class Abc>
struct AlignmentMatrices {
    AlignmentMatrices(const POHmm<Abc>& hmm_x, const POHmm<Abc>& hmm_y)
            : x(hmm_x),
              y(hmm_y),
              F(hmm_x.size() + 1, hmm_y.size() + 1),
              B(hmm_x.size() + 1, hmm_y.size() + 1),
              P(hmm_x.size() + 1, hmm_y.size() + 1),
              S(hmm_x.size() + 1, hmm_y.size() + 1),
              b(hmm_x.size() + 1, hmm_y.size() + 1),
              scale(hmm_x.size() + 1, 1.0),
              ptot(0.0) {}

    POHmm<Abc> x, y;       // PO-HMMs we want to align
    ProbMatrix F;          // forward matrix
    ProbMatrix B;          // backward matrix
    ProbMatrix P;          // posterior probs matrix
    ScoreMatrix S;         // score matrix for MAC algorithm
    BacktraceMatrix b;     // backtrace matrix for MAC algorithm
    Vector<double> scale;  // row scaling factors
    double ptot;           // total forward probability
};


// Algorithm encapsulation for MAC alignment of two PO-HMMs.
template<class Abc>
class PairwiseAligner {
  public:
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::Vertex Vertex;
    typedef typename POHmm<Abc>::OutEdgeIter OutEdgeIter;
    typedef typename POHmm<Abc>::InEdgeIter InEdgeIter;
    typedef typename POHmm<Abc>::VertexIndex VertexIndex;
    typedef typename POHmm<Abc>::VertexVec VertexVec;

    // Constructs an aligner with substitution matrix scoring.
    PairwiseAligner(const SubstitutionMatrix<Abc>* matx, const SubstitutionMatrix<Abc>* maty, double gapf)
            : matx_(matx),
              maty_(maty),
              gapf_(gapf),
              context_score_(0.0)
    { }

    // Constructs an aligner with context scoring. Note that the context library
    // is assumed to be in log-space!
    PairwiseAligner(const SubstitutionMatrix<Abc>* matx,
                    const SubstitutionMatrix<Abc>* maty,
                    const ContextLibrary<Abc>& lib,
                    double gapf,
                    double context_score)
            : matx_(matx),
              maty_(maty),
              priors_(lib.size()),
              cons_(lib.size()),
              gapf_(gapf),
              context_score_(context_score) {
        // Save priors in lin space
        for (size_t k = 0; k < lib.size(); ++k)
            priors_[k] = exp(lib[k].prior);

        // Calculate conservation score of each context
        float bg[] = { 0.3, 0.2, 0.2, 0.3 };  // FIXME: these background freqs work only for DNA!
        for (size_t k = 0; k < lib.size(); ++k) {
            std::vector<double> cons;
            for (size_t i = 0; i < lib.wlen(); ++i) {
                double sum_num = 0.0, sum_den = 0.0;
                for (size_t a = 0; a < Abc::kSize; ++a) {
                    double p = exp(lib[k].probs[i][a]);
                    sum_num += (p * p) / bg[a];
                    sum_den += p / bg[a];
                }
                cons.push_back(log(sum_num) / log(sum_den));
            }
            // Calculate median of conservation values
            std::sort(cons.begin(), cons.end());
            cons_[k] = cons[(lib.wlen() - 1) / 2];
            LOG(ERROR) << strprintf("cons[%s] = %4.2f", lib[k].name.c_str(), cons_[k]);
            //            std::cerr << strprintf("cons[%s] = %4.2f\n", lib[k].name.c_str(), cons_[k]);
        }
    }

    // Constructs an alignment using forward, backward, and MAC algorithm.
    PairAlignment Align(AlignmentMatrices<Abc>& matrices) const;

    // Constructs an alignment using forward, backward, and MAC algorithm.
    PairAlignment Realign(const PairAlignment& ali, AlignmentMatrices<Abc>& matrices) const;

  private:
    // Forward algorithm for two PO-HMMs
    void Forward(AlignmentMatrices<Abc>& matrices) const;

    // Backward algorithm for two PO-HMMs
    void Backward(AlignmentMatrices<Abc>& matrices) const;

    // Calculates MAC alignment on posterior probabilties from forward/backward pass.
    PairAlignment MacAlign(AlignmentMatrices<Abc>& matrices) const;

    // Calculates MAC alignment on posterior probabilties but with the "flawed" normalized recursion.
    PairAlignment MacAlignNormalized(AlignmentMatrices<Abc>& matrices) const;

    // Calculates co-emission probability between context columns 'px' and 'py'
    double ContextProb(const SparseProfileCol& px, const SparseProfileCol& py) const;

    // Calculates co-emission probability between profile columns 'px' and 'py'
    double MatchProb(const ProfileColumn<Abc>& px, const ProfileColumn<Abc>& py) const;

    const SubstitutionMatrix<Abc>* matx_;  // for DNA based scoring at left branch
    const SubstitutionMatrix<Abc>* maty_;  // for DNA based scoring at right branch
    Vector<double> priors_;               // prior probs of contexts (optional)
    Vector<double> cons_;                 // conservation scores of contexts (optional)
    double gapf_;                         // gap factor in MAC algorithm
    double context_score_;                // weight of context score in pairwise match score
};


// Converts a pair state small int into a two character string.
std::string PairStateToString(uint8_t xx);

template<class Abc>
std::ostream& operator<< (std::ostream& out, const AlignmentMatrices<Abc>& m);

std::ostream& operator<< (std::ostream& out, const PairAlignment& ali);

// Extracts matched x-indices in PO-HMM x from pairwise alignment.
inline IndexPath GetPathInX(const PairAlignment& ali) {
    IndexPath rv;
    for (ConstPathIter s = ali.path.begin(); s != ali.path.end(); ++s)
        if (s->pstate == MM || s->pstate == MI)
            rv.push_back(std::make_pair(s->i, s->pstate));
    return rv;
}

// Extracts matched y-indices in PO-HMM y from pairwise alignment.
inline IndexPath GetPathInY(const PairAlignment& ali) {
    IndexPath rv;
    for (ConstPathIter s = ali.path.begin(); s != ali.path.end(); ++s)
        if (s->pstate == MM || s->pstate == IM)
            rv.push_back(std::make_pair(s->j, s->pstate));
    return rv;
}

// Equality operator for steps in an alignment path. Two steps are equal if they
// have the same indices and pair-state. The posterior probability doesn't matter.
inline bool operator== (const Step& lhs, const Step& rhs) {
    return lhs.i == rhs.i && lhs.j == rhs.j && lhs.pstate == rhs.pstate;
}

// Operator 'less' for alignment steps.
inline bool operator< (const Step& lhs, const Step& rhs) {
    if (lhs.i != rhs.i) return lhs.i < rhs.i;
    else if (lhs.j != rhs.j) return lhs.j < rhs.j;
    return lhs.pstate < rhs.pstate;
}

}  // namespace cs

#endif  // CS_PAIRWISE_ALIGNER_H_
