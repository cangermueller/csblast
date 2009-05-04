// Copyright 2009, Andreas Biegert

#ifndef SRC_EMITTER_INL_H_
#define SRC_EMITTER_INL_H_

#include "emitter.h"

#include <cassert>
#include <cmath>

#include "exception.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "log.h"
#include "sequence-inl.h"

namespace cs {

template<class Alphabet>
Emitter<Alphabet>::Emitter(int num_cols, float weight_center, float weight_decay)
    : num_cols_(num_cols),
      center_((num_cols - 1) / 2),
      weights_(1.0f, num_cols) {
  if (num_cols_ % 2 != 1)
    throw Exception("Number of emission columns should be odd but is %i!", num_cols_);
  InitWeights(weight_center, weight_decay);
}

template<class Alphabet>
inline double Emitter<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const CountProfile<Alphabet>& counts,
    int index) const {
  assert(profile.logspace());
  assert(!counts.logspace());

  const int alphabet_size = profile.alphabet_size();
  double rv = 0.0;
  const int beg = std::max(0, index - center_);
  const int end = std::min(counts.num_cols() - 1, index + center_);
  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center_;
    double sum = 0.0;
    for (int a = 0; a < alphabet_size; ++a)
      sum += counts[i][a] * counts.neff(i) * profile[j][a];
    rv += weights_[j] * sum;
  }
  return rv;
}

template<class Alphabet>
inline double Emitter<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const Sequence<Alphabet>& seq,
    int index) const {
  assert(profile.logspace());

  double rv = 0.0;
  const int beg = std::max(0, index - center_);
  const int end = std::min(seq.length() - 1, index + center_);
  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center_;
    rv += weights_[j] * profile[j][seq[i]];
  }
  return rv;
}

template<class Alphabet>
void Emitter<Alphabet>::InitWeights(float weight_center, float weight_decay) {
  weights_[center_] = weight_center;
  for (int i = 1; i <= center_; ++i) {
    float weight = weight_center * pow(weight_decay, i);
    weights_[center_ - i] = weight;
    weights_[center_ + i] = weight;
  }
}

}  // namespace cs

#endif  // SRC_EMITTER_INL_H_
