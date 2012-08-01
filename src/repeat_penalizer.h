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

#ifndef CS_REPEAT_PENALIZER_H_
#define CS_REPEAT_PENALIZER_H_

#include "globals.h"
#include "abstract_state_profile.h"
#include "profile.h"

namespace cs {

// Class implementing CS-BLAST's repeat protein heuristic which penalizes repeat
// regions in a sequence profile by rescaling profile columns.
template<class Alphabet>
class RepeatPenalizer {
 public:
  // Constructs a penalizer object with parameters.
  RepeatPenalizer(float alpha = 0.0f, float beta = 0.1f, float score_min = 8.0f);

  // Rescales repeat regions in given sequence profile.
  void RescaleProfile(const matrix<double>& pp,
                      const float* p_k,
                      Profile<Alphabet>* p) const;

 private:
  // Maximal length of repeats that can be penalized.
  static const int kDelta = 100;
  // Minimal posterior probability for inclusion in abstract state profile.
  static const float kMinPosterior = 0.02f;

  // Baseline penalty to be applied to all profile columns.
  float alpha_;
  // Parameter governing strength of repeat penalty.
  float beta_;
  // Minimum profile-profile score above which repeta penalty takes effect.
  float score_min_;
};  // RepeatPenalizer

}  // namespace cs

#endif  // CS_REPEAT_PENALIZER_H_
