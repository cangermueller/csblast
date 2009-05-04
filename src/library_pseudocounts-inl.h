// Copyright 2009, Andreas Biegert

#ifndef SRC_LIBRARY_PSEUDOCOUNTS_INL_H_
#define SRC_LIBRARY_PSEUDOCOUNTS_INL_H_

#include "library_pseudocounts.h"

#include <cassert>

#include "amino_acid.h"
#include "emitter-inl.h"
#include "count_profile-inl.h"
#include "log.h"
#include "matrix.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils-inl.h"

using std::valarray;

namespace cs {

template<class Alphabet>
inline LibraryPseudocounts<Alphabet>::LibraryPseudocounts(
    const ProfileLibrary<Alphabet>* lib,
    float weight_center,
    float weight_decay)
    : lib_(*lib),
      emitter_(lib->num_cols(), weight_center, weight_decay) {
  assert(lib_.logspace());
}

template<class Alphabet>
void LibraryPseudocounts<Alphabet>::add_to_sequence(
    const Sequence<Alphabet>& seq,
    const Admixture& pca,
    Profile<Alphabet>* profile) const {
  LOG(DEBUG1) << "Adding context-specific library pseudocounts to sequence ...";
  LOG(DEBUG1) << seq;

  if (seq.length() != profile->num_cols())
    throw Exception("Cannot add context-specific pseudocounts: "
                    "sequence and profile have different length!");

  const int length        = seq.length();
  const int num_profiles  = lib_.num_profiles();
  const int alphabet_size = Alphabet::instance().size();
  const int center        = lib_.center();
  const float tau         = pca(1.0f);  // number of effective seqs is one

  // Profile probabilities P(p_k|X_i) at position i
  valarray<double> prob(0.0, num_profiles);
  // Pseudocount vector P(a|X_i) at position i
  valarray<double> pc(0.0, alphabet_size);
  // Output profile with pseudocounts
  Profile<Alphabet>& p = *profile;

  for (int i = 0; i < length; ++i) {
    // Calculate profile probabilities P(p_k|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      prob[k] = fast_pow2(emitter_(lib_[k], seq, i)) * lib_[k].prior();
    }
    prob /= prob.sum();  // normalization

    // Calculate pseudocount vector P(a|X_i)
    pc = 0.0;
    for (int k = 0; k < num_profiles; ++k) {
      for(int a = 0; a < alphabet_size; ++a) {
        pc[a] += prob[k] * fast_pow2(lib_[k][center][a]);
      }
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence
    for(int a = 0; a < alphabet_size; ++a) {
      assert(pc[a] > 0.0);  // should be always true because of log-scaling
      const float pa = (1.0f - tau) *
        (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * pc[a];
      p[i][a] = p.logspace() ? fast_log2(pa) : pa;
    }
  }
  normalize(profile);

  LOG(DEBUG1) << *profile;
}

template<class Alphabet>
void LibraryPseudocounts<Alphabet>::add_to_profile(
    const Admixture& pca,
    CountProfile<Alphabet>* profile) const {
  assert(!profile->logspace());

  LOG(DEBUG2) << "Adding context-specific library pseudocounts to profile ...";
  LOG(DEBUG2) << *profile;

  const int length        = profile->num_cols();
  const int num_profiles  = lib_.num_profiles();
  const int center        = lib_.center();
  const int alphabet_size = Alphabet::instance().size();

  // Unscaled log emission probabilities P(X_i|p_k)
  valarray<double> log_ep(0.0, num_profiles);
  // Profile probabilities P(p_k|X_i) at position i
  valarray<double> prob(0.0, num_profiles);
  // Pseudocount matrix P(a|X_i)
  matrix<double> pc(length, alphabet_size, 0.0);
  // Output profile with pseudocounts
  CountProfile<Alphabet>& p = *profile;

  for (int i = 0; i < length; ++i) {
    // Calculate unscaled emission probs
    for (int k = 0; k < num_profiles; ++k) {
      log_ep[k] = emitter_(lib_[k], p, i);
    }
    double log_ep_average = log_ep.sum() / num_profiles;
    log_ep -= log_ep_average;  // scale logs

    // Calculate profile probabilities P(p_k|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      prob[k] = fast_pow2(log_ep[k]) * lib_[k].prior();
    }
    prob /= prob.sum();  // normalization

    // Calculate pseudocount vector P(a|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      for(int a = 0; a < alphabet_size; ++a) {
        pc[i][a] += prob[k] * fast_pow2(lib_[k][center][a]);
      }
    }
    normalize_to_one(&pc[i][0], alphabet_size);
  }

  // Add pseudocounts to profile
  for (int i = 0; i < length; ++i) {
    const float tau = pca(p.neff(i));

    for(int a = 0; a < alphabet_size; ++a) {
      assert(pc[i][a] > 0.0);  // should be always true because of log-scaling
      p[i][a] = (1.0f - tau) * p[i][a] + tau * pc[i][a];
    }
  }
  normalize(profile);

  LOG(DEBUG2) << *profile;
}

}  // namespace cs

#endif  // SRC_LIBRARY_PSEUDOCOUNTS_INL_H_
