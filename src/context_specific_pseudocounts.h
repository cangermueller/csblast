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

#ifndef CS_CONTEXT_SPECIFIC_PSEUDOCOUNTS_H_
#define CS_CONTEXT_SPECIFIC_PSEUDOCOUNTS_H_

#include <stdio.h>
#include <stdlib.h>

#include "count_profile-inl.h"
#include "matrix.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
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
