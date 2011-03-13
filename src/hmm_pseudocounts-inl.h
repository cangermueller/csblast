// Copyright 2009, Andreas Biegert

#ifndef CS_HMM_PSEUDOCOUNTS_INL_H_
#define CS_HMM_PSEUDOCOUNTS_INL_H_

#include "hmm_pseudocounts.h"

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "count_profile-inl.h"
#include "mult_emission-inl.h"
#include "forward_backward_algorithm-inl.h"
#include "hmm-inl.h"
#include "log.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils.h"

namespace cs {

template<class Alphabet>
HMMPseudocounts<Alphabet>::HMMPseudocounts(const Hmm<Alphabet>* hmm,
                                           int window_length,
                                           float weight_center,
                                           float weight_decay)
    : hmm_(*hmm),
      emission_(window_length, weight_center, weight_decay) {
  assert(hmm_.states_logspace());
  assert(!hmm_.transitions_logspace());
}

template<class Alphabet>
void HMMPseudocounts<Alphabet>::AddPseudocountsToSequence(
    const Sequence<Alphabet>& seq,
    const Admixture& pca,
    Profile<Alphabet>* profile) const {
  LOG(DEBUG2) << "Adding context-specific HMM pseudocounts to sequence ...";
  LOG(DEBUG2) << seq;

  if (seq.length() != profile->num_cols())
    throw Exception("Cannot add context-specific pseudocounts: "
                    "sequence and profile have different length!");

  const int length        = seq.length();
  const int num_states    = hmm_.num_states();
  const int alphabet_size = Alphabet::instance().size();
  const int center        = hmm_.center();
  const float tau         = pca(1.0f);  // number of effective seqs is one

  // Pseudocount vector P(a|X_i) at position i
  std::valarray<float> pc(0.0f, alphabet_size);
  // Output profile with pseudocounts
  Profile<Alphabet>& p = *profile;

  // Calculate posterior state probabilities with forward-backward algorithm
  ForwardBackwardMatrices fbm(length, hmm_.num_states());
  ForwardBackwardAlgorithm(hmm_, seq, emission_, &fbm);

  for (int i = 0; i < length; ++i) {
    // Calculate pseudocount vector P(a|X_i)
    for(int a = 0; a < alphabet_size; ++a) {
      pc[a] = 0.0f;
      for (int k = 0; k < num_states; ++k) {
        pc[a] += fbm.f[i][k] * fbm.b[i][k] * fast_pow2(hmm_[k][center][a]);
      }
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence by storing probabilities in output profile
    for(int a = 0; a < alphabet_size; ++a) {
      float pa = tau * pc[a] + (1.0f - tau) *
        (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f);
      p[i][a] = p.logspace() ? fast_log2(pa) : pa;
    }
  }

  Normalize(profile);
  RescaleProfile(fbm.pp, profile);

  LOG(DEBUG2) << *profile;
}

template<class Alphabet>
void HMMPseudocounts<Alphabet>::AddPseudocountsToProfile(
    const Admixture& pca,
    CountProfile<Alphabet>* profile) const {
  assert(!profile->logspace());

  LOG(DEBUG2) << "Adding context-specific, HMM-derived pseudocounts to profile ...";
  LOG(DEBUG2) << *profile;

  const int length        = profile->num_cols();
  const int num_states    = hmm_.num_states();
  const int center        = hmm_.center();
  const int alphabet_size = Alphabet::instance().size();

  // Pseudocount vector P(a|X_i) at position i
  std::valarray<float> pc(0.0f, alphabet_size);
  // Output profile with pseudocounts
  CountProfile<Alphabet>& p = *profile;

  // Calculate posterior state probabilities with forward-backward algorithm
  ForwardBackwardMatrices fbm(length, hmm_.num_states());
  ForwardBackwardAlgorithm(hmm_, p, emission_, &fbm);

  for (int i = 0; i < length; ++i) {
    // Calculate pseudocount vector P(a|X_i)
    for(int a = 0; a < alphabet_size; ++a) {
      pc[a] = 0.0f;
      for (int k = 0; k < num_states; ++k) {
        pc[a] += fbm.pp[i][k] * fast_pow2(hmm_[k][center][a]);
      }
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence by storing probabilities in output profile
    float tau = pca(p.neff(i));
    for(int a = 0; a < alphabet_size; ++a) {
      p[i][a] = (1.0f - tau) * p[i][a] + tau * pc[a];
    }
  }

  Normalize(profile);
  RescaleProfile(fbm.pp, profile);

  LOG(DEBUG2) << *profile;
}

template<class Alphabet>
void HMMPseudocounts<Alphabet>::CalculatePosteriorsForSequence(
    const Sequence<Alphabet>& seq,
    matrix<double>* m) const {
  const int length        = seq.length();
  const int num_states    = hmm_.num_states();

  // Calculate posterior state probabilities with forward-backward algorithm
  ForwardBackwardMatrices fbm(length, hmm_.num_states());
  ForwardBackwardAlgorithm(hmm_, seq, emission_, &fbm);

  matrix<double>& pp = *m;
  for (int i = 0; i < length; ++i)
    for (int k = 0; k < num_states; ++k)
      pp[i][k] = fbm.pp[i][k];
}

template<class Alphabet>
void HMMPseudocounts<Alphabet>::CalculatePosteriorsForProfile(
    const CountProfile<Alphabet>& profile,
    matrix<double>* m) const {
  const int length        = profile.num_cols();
  const int num_states    = hmm_.num_states();

  // Calculate posterior state probabilities with forward-backward algorithm
  ForwardBackwardMatrices fbm(length, hmm_.num_states());
  ForwardBackwardAlgorithm(hmm_, profile, emission_, &fbm);

  matrix<double>& pp = *m;
  for (int i = 0; i < length; ++i)
    for (int k = 0; k < num_states; ++k)
      pp[i][k] = fbm.pp[i][k];
}

template<class Alphabet>
void HMMPseudocounts<Alphabet>::RescaleProfile(const matrix<double>& pp,
                                               Profile<Alphabet>* profile) const {
  if (penalizer_) {
    valarray<float> p_k(0.0, pp.num_cols());
    for (int k = 0; k < pp.num_cols(); ++k)
      for (int i = 0; i < pp.num_rows(); ++i)
        p_k[k] += pp[i][k];
    penalizer_->RescaleProfile(pp, &p_k[0], profile);
  }
}

}  // namespace cs

#endif  // CS_HMM_PSEUDOCOUNTS_INL_H_
