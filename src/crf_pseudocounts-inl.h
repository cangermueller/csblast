// Copyright 2009, Andreas Biegert

#ifndef CS_CRF_PSEUDOCOUNTS_INL_H_
#define CS_CRF_PSEUDOCOUNTS_INL_H_

#include "crf_pseudocounts.h"

namespace cs {

template<class Abc>
CrfPseudocounts<Abc>::CrfPseudocounts(const Crf<Abc>& crf) : crf_(crf) {}

template<class Abc>
void CrfPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq,
                                         const Admix& pca,
                                         Profile<Abc>& p) const {
  assert_eq(seq.length(), p.length());
  LOG(INFO) << "Adding CRF pseudocounts to sequence ...";

  const size_t center = crf_.center();
  const double tau = pca(1.0);  // effective number of sequences is one
  Matrix<double> pp(seq.length(), crf_.size(), 0.0);  // posterior probabilities
  Vector<double> pc(Abc::kSize);                      // pseudocount vector P(a|X_i)

  // Calculate and add pseudocounts for each sequence window X_i separately
  for (size_t i = 0; i < seq.length(); ++i) {
    double* ppi = &pp[i][0];

    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    double max = -DBL_MAX;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = crf_[k].bias_weight +
        ContextScore(crf_[k].context_weights, seq, i, center);
      if (ppi[k] > max)
        max = ppi[k];  // needed for log-sum-exp trick
    }

    // Log-sum-exp trick begins here
    double sum = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k)
      sum += exp(ppi[k] - max);
    double tmp = max + log(sum);
    Assign(pc, 0.0);
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = exp(ppi[k] - tmp);
      // Calculate pseudocount vector P(a|X_i)
      for(size_t a = 0; a < Abc::kSize; ++a)
        pc[a] += ppi[k] * crf_[k].pc[a];
    }
    Normalize(&pc[0], Abc::kSize);  // FIXME: is this really needed?

    // Add pseudocounts to sequence
    for(size_t a = 0; a < Abc::kSize; ++a)
      p[i][a] = (1.0 - tau) * (seq[i] == a ? 1.0 : 0.0) + tau * pc[a];
  }
}

template<class Abc>
void CrfPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp,
                                        const Admix& pca,
                                        Profile<Abc>& p) const {
  assert_eq(cp.counts.length(), p.length());
  LOG(INFO) << "Adding library pseudocounts to profile ...";

  const size_t center = crf_.center();
  Matrix<double> pp(cp.counts.length(), crf_.size(), 0.0);  // posterior probs
  Vector<double> pc(Abc::kSize);  // pseudocount vector P(a|X_i)

  // Calculate and add pseudocounts for each sequence window X_i separately
  for (size_t i = 0; i < cp.counts.length(); ++i) {
    double* ppi = &pp[i][0];

    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    double max = -DBL_MAX;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = crf_[k].bias_weight +
        ContextScore(crf_[k].context_weights, cp, i, center);
      if (ppi[k] > max)
        max = ppi[k];  // needed for log-sum-exp trick
    }

     // Log-sum-exp trick begins here
    double sum = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k)
      sum += exp(ppi[k] - max);
    double tmp = max + log(sum);
    Assign(pc, 0.0);
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = exp(ppi[k] - tmp);
      // Calculate pseudocount vector P(a|X_i)
      for(size_t a = 0; a < Abc::kSize; ++a)
        pc[a] += ppi[k] * crf_[k].pc[a];
    }
    Normalize(&pc[0], Abc::kSize);  // FIXME: is this really needed?

    // Add pseudocounts to profile
    for(size_t a = 0; a < Abc::kSize; ++a) {
      double tau = pca(cp.neff[i]);
      p[i][a] = (1.0 - tau) * cp.counts[i][a] / cp.neff[i] + tau * pc[a];
    }
  }
}

}  // namespace cs

#endif  // CS_LIBRARY_PSEUDOCOUNTS_INL_H_
