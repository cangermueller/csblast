// Copyright 2009, Andreas Biegert

#ifndef SRC_EMITTER_INL_H_
#define SRC_EMITTER_INL_H_

#include "mult_emission.h"

#include <cassert>
#include <cmath>

#include "amino_acid.h"
#include "exception.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "log.h"
#include "sequence-inl.h"

namespace cs {

template<class Alphabet>
MultEmission<Alphabet>::MultEmission(int num_cols,
                                     float weight_center,
                                     float weight_decay)
    : num_cols_(num_cols),
      center_((num_cols - 1) / 2),
      weights_(1.0f, num_cols) {
  if (num_cols_ % 2 != 1)
    throw Exception("Number of emission columns should be odd but is %i!", num_cols_);
  InitWeights(weight_center, weight_decay);
}

template<class Alphabet>
inline double MultEmission<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const CountProfile<Alphabet>& counts,
    int index) const {
  assert(profile.logspace());
  assert(!counts.logspace());

  const int alphabet_size = profile.alphabet_size();
  const int beg = std::max(0, index - center_);
  const int end = std::min(counts.num_cols() - 1, index + center_);
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center_;
    double sum = 0.0;
    for (int a = 0; a < alphabet_size; ++a)
      sum += counts[i][a] * counts.neff(i) * profile[j][a];
    rv += weights_[j] * sum;
  }

  return rv;
}

template<>
inline double MultEmission<AminoAcid>::operator() (
    const ContextProfile<AminoAcid>& profile,
    const CountProfile<AminoAcid>& counts,
    int index) const {
  assert(profile.logspace());
  assert(!counts.logspace());

  const int beg = std::max(0, index - center_);
  const int end = std::min(counts.num_cols() - 1, index + center_);
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center_;
    double sum = 0.0;
    sum += counts[i][0] * counts.neff(i) * profile[j][0];
    sum += counts[i][1] * counts.neff(i) * profile[j][1];
    sum += counts[i][2] * counts.neff(i) * profile[j][2];
    sum += counts[i][3] * counts.neff(i) * profile[j][3];
    sum += counts[i][4] * counts.neff(i) * profile[j][4];
    sum += counts[i][5] * counts.neff(i) * profile[j][5];
    sum += counts[i][6] * counts.neff(i) * profile[j][6];
    sum += counts[i][7] * counts.neff(i) * profile[j][7];
    sum += counts[i][8] * counts.neff(i) * profile[j][8];
    sum += counts[i][9] * counts.neff(i) * profile[j][9];
    sum += counts[i][10] * counts.neff(i) * profile[j][10];
    sum += counts[i][11] * counts.neff(i) * profile[j][11];
    sum += counts[i][12] * counts.neff(i) * profile[j][12];
    sum += counts[i][13] * counts.neff(i) * profile[j][13];
    sum += counts[i][14] * counts.neff(i) * profile[j][14];
    sum += counts[i][15] * counts.neff(i) * profile[j][15];
    sum += counts[i][16] * counts.neff(i) * profile[j][16];
    sum += counts[i][17] * counts.neff(i) * profile[j][17];
    sum += counts[i][18] * counts.neff(i) * profile[j][18];
    sum += counts[i][19] * counts.neff(i) * profile[j][19];

    rv += weights_[j] * sum;
  }

  return rv;
}

template<class Alphabet>
inline double MultEmission<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const Sequence<Alphabet>& seq,
    int index) const {
  assert(profile.logspace());

  const int beg = std::max(0, index - center_);
  const int end = std::min(seq.length() - 1, index + center_);
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center_;
    rv += weights_[j] * profile[j][seq[i]];
  }

  return rv;
}

template<class Alphabet>
void MultEmission<Alphabet>::InitWeights(float weight_center, float weight_decay) {
  weights_[center_] = weight_center;
  for (int i = 1; i <= center_; ++i) {
    float weight = weight_center * pow(weight_decay, i);
    weights_[center_ - i] = weight;
    weights_[center_ + i] = weight;
  }
}

}  // namespace cs

#endif  // SRC_EMITTER_INL_H_
