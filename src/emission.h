// Copyright 2009, Andreas Biegert

#ifndef CS_EMISSION_H_
#define CS_EMISSION_H_

#include "count_profile-inl.h"
#include "po_hmm-inl.h"
#include "profile-inl.h"
#include "profile_column.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Functor for calculating multinomial emission probabilities for context profiles.
template<class Abc>
class Emission {
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::Vertex Vertex;
    typedef typename POHmm<Abc>::OutEdgeIter OutEdgeIter;
    typedef typename POHmm<Abc>::InEdgeIter InEdgeIter;
    typedef typename POHmm<Abc>::VertexIndex VertexIndex;

  public:

    // Constructs an emission object with positional window weights.
    Emission(size_t wlen,
	     double w_center = 1.6,
	     double w_decay = 0.85,
	     const SubstitutionMatrix<Abc>* sm = NULL)
	    : center_((wlen - 1) / 2), weights_(wlen), logp_(0.0) {
	assert(wlen & 1);

	weights_[center_] = w_center;
	for (size_t i = 1; i <= center_; ++i) {
	    double weight = w_center * pow(w_decay, i);
	    weights_[center_ - i] = weight;
	    weights_[center_ + i] = weight;
	}

	if (sm) {
	    for (size_t a = 0; a < Abc::kSize; ++a)
		logp_[a] = log(sm->p(a));
	}
    }

    // Calculates the sum of positional weights.
    float GetSumWeights() const {
	double sum = 0.0;
	for (size_t i = 0; i < weights_.size(); ++i) sum += weights_[i];
	return sum;
    }

    // Calculates the log of the probability that profile 'p' emits the sequence
    // window centered at index 'idx' in 'seq'. Note that the normalization factor
    // that is usualy used in multinomial distributions is left out since it
    // cancels out anyway.
    double operator() (const Profile<Abc>& p,
		       const Sequence<Abc>& seq,
		       size_t idx) const {
	assert(p.length() & 1);
	assert_eq(weights_.size(), p.length());

	const size_t beg = MAX(0, static_cast<int>(idx - center_));
	const size_t end = MIN(seq.length(), idx + center_ + 1);
	double rv = 0.0;
	for(size_t i = beg, j = beg - idx + center_; i < end; ++i, ++j) {
	    rv += weights_[j] * (p[j][seq[i]] - logp_[seq[i]]);
	}
	return rv;
    }

    // Calculates the log of the probability that profile 'p' emits the counts in
    // a window centered at index 'idx' in 'c'. Note that the normalization factor
    // that is usualy used in multinomial distributions is left out since it
    // cancels out anyway.
    double operator() (const Profile<Abc>& p,
		       const CountProfile<Abc>& cp,
		       size_t idx) const {
	assert(p.length() & 1);
	assert_eq(weights_.size(), p.length());

	const size_t beg = MAX(0, static_cast<int>(idx - center_));
	const size_t end = MIN(cp.counts.length(), idx + center_ + 1);
	double sum, rv = 0.0;
	for(size_t i = beg, j = beg - idx + center_; i < end; ++i, ++j) {
	    sum = 0.0;
	    for (size_t a = 0; a < Abc::kSize; ++a)
		sum += cp.counts[i][a] * (p[j][a] - logp_[a]);
	    rv += weights_[j] * sum;
	}
	return rv;
    }

    // Calculates the log of the probability that profile 'p' emits the counts in
    // a window centered around POG vertex 'v'. Note that the normalization factor
    // that is usualy used in multinomial distributions is left out since it
    // cancels out anyway.
    double operator() (const Profile<Abc>& p,
		       const Graph& g,
		       Vertex v) const {
	assert(p.length() & 1);
	assert_eq(weights_.size(), p.length());
	VertexIndex index;
	// Calculate emission at central vertex 'v'
	double rv = 0.0, sum = 0.0;
	for (size_t a = 0; a < Abc::kSize; ++a)
	    sum += g[v].counts[a] * (p[center_][a] - logp_[a]);
	rv += weights_[center_] * sum;
	// Collect emission probs within context window in forwad direction
	if (center_ != 0) {
	    OutEdgeIter ei, edge_end;
	    for (tie(ei, edge_end) = out_edges(v, g); ei != edge_end; ++ei) {
		if (index[target(*ei, g)] != kStartEndVertex)
		    rv += g[*ei].weight * ScoreFwd(p, g, target(*ei, g), center_ + 1);
	    }
	}
	// Collect emission probs within context window in backward direction
	if (center_ != 0) {
	    InEdgeIter ei, edge_end;
	    for (tie(ei, edge_end) = in_edges(v, g); ei != edge_end; ++ei) {
		if (index[source(*ei, g)] != kStartEndVertex)
		    rv += g[*ei].weight_rev * ScoreBwd(p, g, source(*ei, g), center_ - 1);
	    }
	}
	return rv;
    }

  private:
    double ScoreFwd(const Profile<Abc>& p,
		    const Graph& g,
		    Vertex v,
		    int j) const {
	assert(j < static_cast<int>(weights_.size()));
	VertexIndex index;
	// Calculate emission score for counts at vertex 'v'
	double rv = 0.0, sum = 0.0;
	for (size_t a = 0; a < Abc::kSize; ++a)
	    sum += g[v].counts[a] * (p[j][a] - logp_[a]);
	rv += weights_[j] * sum;
	// Collect emission probs within context window in forwad direction
	if (j + 1 < static_cast<int>(weights_.size())) {
	    OutEdgeIter ei, edge_end;
	    for (tie(ei, edge_end) = out_edges(v, g); ei != edge_end; ++ei) {
		if (index[target(*ei, g)] != kStartEndVertex)
		    rv += g[*ei].weight * ScoreFwd(p, g, target(*ei, g), j + 1);
	    }
	}
	return rv;
    }

    double ScoreBwd(const Profile<Abc>& p,
		    const Graph& g,
		    Vertex v,
		    int j) const {
	assert(j >= 0);
	VertexIndex index;
	// Calculate emission score for counts at vertex 'v'
	double rv = 0.0, sum = 0.0;
	for (size_t a = 0; a < Abc::kSize; ++a)
	    sum += g[v].counts[a] * (p[j][a] - logp_[a]);
	rv += weights_[j] * sum;
	// Collect emission probs within context window in backward direction
	if (j - 1 >= 0) {
	    InEdgeIter ei, edge_end;
	    for (tie(ei, edge_end) = in_edges(v, g); ei != edge_end; ++ei) {
		if (index[source(*ei, g)] != kStartEndVertex)
		    rv += g[*ei].weight_rev * ScoreBwd(p, g, source(*ei, g), j - 1);
	    }
	}
	return rv;
    }

    size_t center_;            // index of central column in context window
    Vector<double> weights_;   // positional window weights
    ProfileColumn<Abc> logp_;  // log of background frequencies
};  // class Emission

}  // namespace cs

#endif  // CS_EMISSION_H_
