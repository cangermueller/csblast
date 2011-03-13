// Copyright 2009, Andreas Biegert

#ifndef CS_CONTEXT_SPECIFIC_PSEUDOCOUNTS_H_
#define CS_CONTEXT_SPECIFIC_PSEUDOCOUNTS_H_

#include <stdio.h>
#include <stdlib.h>

#include "count_profile-inl.h"
#include "matrix.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "repeat_penalizer.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from a library
// of context profiles.
template<class Alphabet>
class ContextSpecificPseudocounts : public Pseudocounts<Alphabet> {
  public:
    ContextSpecificPseudocounts() : penalizer_(NULL) {}
    virtual ~ContextSpecificPseudocounts() {}

    // Sets repeat penalizer
    void set_repeat_penalizer(const RepeatPenalizer<Alphabet>* penalizer) {
	penalizer_ = penalizer;
    }
    // Rescales sequence profile with repeat penalty heuristic
    virtual void RescaleProfile(const matrix<double>& pp,
				Profile<Alphabet>* profile) const = 0;
    // Calculates posterior probability matrix for an input sequence.
    virtual void CalculatePosteriorsForSequence(const Sequence<Alphabet>& seq, matrix<double>* m) const = 0;

    // Calculates posterior probability matrix for an input profile.
    virtual void CalculatePosteriorsForProfile(const CountProfile<Alphabet>& profile, matrix<double>* m) const = 0;

  protected:
    // Pointer to repeat penalizer for rescaling output profile in CS-BLAST
    const RepeatPenalizer<Alphabet>* penalizer_;

  private:
    DISALLOW_COPY_AND_ASSIGN(ContextSpecificPseudocounts);
};  // ContextSpecificPseudocounts

}  // namespace cs

#endif  // CS_CONTEXT_SPECIFIC_PSEUDOCOUNTS_H_
