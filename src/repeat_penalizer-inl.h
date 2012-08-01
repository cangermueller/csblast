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

#ifndef CS_REPEAT_PENALIZER_INL_H_
#define CS_REPEAT_PENALIZER_INL_H_

#include "repeat_penalizer.h"

#include "globals.h"
#include "abstract_state_profile.h"
#include "log.h"
#include "matrix.h"
#include "profile-inl.h"

namespace cs {

template<class Alphabet>
RepeatPenalizer<Alphabet>::RepeatPenalizer(float alpha,
                                           float beta,
                                           float score_min)
    : alpha_(alpha),
      beta_(beta),
      score_min_(score_min) {}

template<class Alphabet>
void RepeatPenalizer<Alphabet>::RescaleProfile(const matrix<double>& pp,
                                               const float* p_k,
                                               Profile<Alphabet>* p) const {
  typedef AbstractStateProfile::ColIter ColIter;
  typedef AbstractStateProfile::Element Element;

  const int length = pp.num_rows();
  const int alphabet_size = Alphabet::instance().size();
  AbstractStateProfile asp(pp, kMinPosterior);
  Profile<Alphabet>& profile = *p;

  for (int i = 0; i < length; ++i) {
    const int beg = MAX(i - kDelta, 0);
    const int end = MIN(i + kDelta, length - 1);
    float sum_delta = 0.0f;

    for(int j = beg; j <= end; ++j) {
      if (j != i) {
        float s = 0.0f;
        // Calculate profile-profile score of column i with column (i+delta).
        for (ColIter it = asp.col_begin(i); it != asp.col_end(i); ++it) {
          s += it->prob * asp[j].get(it->index).prob / p_k[it->index];
        }
        if (s > 0.0f) {
          s = fast_log2(s);
          if (s >= score_min_)
            sum_delta += s;
        }
      }
    }

    float scale_penalty = fast_pow2(-(alpha_ + beta_ * sqrt(sum_delta)));
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] *= scale_penalty;
    LOG(INFO) << strprintf("i=%-4i  scale=%.2f", i, scale_penalty);
  }
}

}  // namespace cs

#endif  // CS_REPEAT_PENALIZER_INL_H_
