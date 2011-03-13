// Copyright 2009, Andreas Biegert

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
