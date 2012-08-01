/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_MULT_EMISSION_INL_H_
#define CS_MULT_EMISSION_INL_H_

#include "mult_emission.h"

#include <math.h>

#include <algorithm>

#include "globals.h"
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
    throw Exception("Number of emission columns should be odd but is %i!",
                    num_cols_);
  InitWeights(weight_center, weight_decay);
}

template<class Alphabet>
inline double MultEmission<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const CountProfile<Alphabet>& count_profile,
    int index) const {
  assert(profile.logspace());
  assert(!count_profile.logspace());

  const int alphabet_size = profile.alphabet_size();
  const int beg = MAX(0, index - center_);
  const int end = MIN(count_profile.num_cols() - 1, index + center_);
  const int prof_center = profile.center();
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + prof_center;
    double sum = 0.0;
    for (int a = 0; a < alphabet_size; ++a)
      sum += count_profile.counts(i, a) * profile[j][a];
    rv += sum * weights_[i - index + center_];
  }

  return rv;
}

template<>
inline double MultEmission<AminoAcid>::operator() (
    const ContextProfile<AminoAcid>& profile,
    const CountProfile<AminoAcid>& count_profile,
    int index) const {
  assert(profile.logspace());
  assert(!count_profile.logspace());

  const int beg = MAX(0, index - center_);
  const int end = MIN(count_profile.num_cols() - 1, index + center_);
  const int prof_center = profile.center();
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + prof_center;
    double sum = 0.0;
    sum += count_profile.counts(i,0) * profile[j][0];
    sum += count_profile.counts(i,1) * profile[j][1];
    sum += count_profile.counts(i,2) * profile[j][2];
    sum += count_profile.counts(i,3) * profile[j][3];
    sum += count_profile.counts(i,4) * profile[j][4];
    sum += count_profile.counts(i,5) * profile[j][5];
    sum += count_profile.counts(i,6) * profile[j][6];
    sum += count_profile.counts(i,7) * profile[j][7];
    sum += count_profile.counts(i,8) * profile[j][8];
    sum += count_profile.counts(i,9) * profile[j][9];
    sum += count_profile.counts(i,10) * profile[j][10];
    sum += count_profile.counts(i,11) * profile[j][11];
    sum += count_profile.counts(i,12) * profile[j][12];
    sum += count_profile.counts(i,13) * profile[j][13];
    sum += count_profile.counts(i,14) * profile[j][14];
    sum += count_profile.counts(i,15) * profile[j][15];
    sum += count_profile.counts(i,16) * profile[j][16];
    sum += count_profile.counts(i,17) * profile[j][17];
    sum += count_profile.counts(i,18) * profile[j][18];
    sum += count_profile.counts(i,19) * profile[j][19];

    rv += sum * weights_[i - index + center_];
  }

  return rv;
}

template<class Alphabet>
inline double MultEmission<Alphabet>::operator() (
    const ContextProfile<Alphabet>& profile,
    const Sequence<Alphabet>& seq,
    int index) const {
  assert(profile.logspace());

  const int beg = MAX(0, index - center_);
  const int end = MIN(seq.length() - 1, index + center_);
  const int prof_center = profile.center();
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + prof_center;
    rv += weights_[i - index + center_] * profile[j][seq[i]];
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

#endif  // CS_MULT_EMISSION_INL_H_
