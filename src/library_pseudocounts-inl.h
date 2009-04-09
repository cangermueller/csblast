// Copyright 2009, Andreas Biegert

#ifndef SRC_LIBRARY_PSEUDOCOUNTS_INL_H_
#define SRC_LIBRARY_PSEUDOCOUNTS_INL_H_

#include "library_pseudocounts.h"

#include <cassert>
#include <cmath>

#include "count_profile-inl.h"
#include "log.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils.h"

namespace cs {

template<class Alphabet>
inline LibraryPseudocounts<Alphabet>::LibraryPseudocounts(
    const ProfileLibrary<Alphabet>* lib,
    const EmissionParams& params)
    : lib_(*lib),
      emitter_(lib->num_cols(), params) {
  assert(lib_.logspace());
}

template<class Alphabet>
void LibraryPseudocounts<Alphabet>::add_to_sequence(
    const Sequence<Alphabet>& seq,
    const Admixture& pca,
    Profile<Alphabet>* profile) const {
  LOG(DEBUG2) << "Adding context-specific library pseudocounts to sequence ...";
  LOG(DEBUG2) << *profile;

  if (seq.length() != profile->num_cols())
    throw Exception("Cannot add context-specific pseudocounts: "
                    "sequence and profile have different length!");

  const int length        = seq.length();
  const int num_profiles  = lib_.num_profiles();
  const int alphabet_size = Alphabet::instance().size();
  const int center        = lib_.center();
  const float tau         = pca(1.0f);  // number of effective seqs is one

  // Profile probabilities P(p_k|X_i) at position i
  std::valarray<float> prob(0.0f, num_profiles);
  // Pseudocount vector P(a|X_i) at position i
  std::valarray<float> pc(0.0f, alphabet_size);
  // Output profile with pseudocounts
  Profile<Alphabet>& p = *profile;

  for (int i = 0; i < length; ++i) {
    // Calculate profile probabilities P(p_k|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      prob[k] = lib_[k].prior() * pow(2.0, emitter_(lib_[k], seq, i));
    }
    prob /= prob.sum();  // normalization

    // Calculate pseudocount vector P(a|X_i)
    for(int a = 0; a < alphabet_size; ++a) {
      pc[a] = 0.0f;
      for (int k = 0; k < num_profiles; ++k) {
        pc[a] += prob[k] * pow(2.0, lib_[k][center][a]);
      }
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence by storing probabilities in output profile
    for(int a = 0; a < alphabet_size; ++a) {
      float pa = (1.0f - tau) * (static_cast<int>(seq[i]) == a ?
                                 1.0f : 0.0f) + tau * pc[a];
      p[i][a] = p.logspace() ? log2(pa) : pa;
    }
  }
  normalize(profile);

  LOG(DEBUG2) << *profile;
}

template<class Alphabet>
void LibraryPseudocounts<Alphabet>::add_to_profile(
    const Admixture& pca,
    CountProfile<Alphabet>* profile) const {
  assert(!profile->has_counts());
  assert(!profile->logspace());

  LOG(DEBUG2) << "Adding context-specific library pseudocounts to profile ...";
  LOG(DEBUG2) << *profile;

  const int length        = profile->num_cols();
  const int num_profiles  = lib_.num_profiles();
  const int center        = lib_.center();
  const int alphabet_size = Alphabet::instance().size();

  // Profile probabilities P(p_k|X_i) at position i
  std::valarray<float> prob(0.0f, num_profiles);
  // Pseudocount vector P(a|X_i) at position i
  std::valarray<float> pc(0.0f, alphabet_size);
  // Output profile with pseudocounts
  CountProfile<Alphabet>& p = *profile;

  for (int i = 0; i < length; ++i) {
    // Calculate profile probabilities P(p_k|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      prob[k] = lib_[k].prior() * pow(2.0, emitter_(lib_[k], p, i));
    }
    prob /= prob.sum();  // normalization

    // Calculate pseudocount vector P(a|X_i)
    for(int a = 0; a < alphabet_size; ++a) {
      pc[a] = 0.0f;
      for (int k = 0; k < num_profiles; ++k) {
        pc[a] += prob[k] * pow(2.0, lib_[k][center][a]);
      }
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence by storing probabilities in output profile
    float tau = pca(p.neff(i));
    for(int a = 0; a < alphabet_size; ++a) {
      p[i][a] = (1.0f - tau) * p[i][a] + tau * pc[a];
    }
  }
  normalize(profile);

  LOG(DEBUG2) << *profile;
}

}  // namespace cs

#endif  // SRC_LIBRARY_PSEUDOCOUNTS_INL_H_
