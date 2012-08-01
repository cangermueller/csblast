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

#ifndef CS_HMM_PSEUDOCOUNTS_H_
#define CS_HMM_PSEUDOCOUNTS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <valarray>

#include "count_profile-inl.h"
#include "mult_emission.h"
#include "hmm-inl.h"
#include "matrix.h"
#include "profile-inl.h"
#include "context_specific_pseudocounts.h"
#include "repeat_penalizer.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from context HMM.
template<class Alphabet>
class HMMPseudocounts : public ContextSpecificPseudocounts<Alphabet> {
 public:
  HMMPseudocounts(const Hmm<Alphabet>* hmm,
                  int window_length,
                  float weight_center,
                  float weight_decay);
  virtual ~HMMPseudocounts() {}

  // Adds context-specific pseudocounts to sequence and stores resulting
  // frequencies in given profile.
  virtual void AddPseudocountsToSequence(const Sequence<Alphabet>& seq,
                               const Admixture& pca,
                               Profile<Alphabet>* profile) const;
  // Adds context-specific pseudocounts to alignment derived profile.
  virtual void AddPseudocountsToProfile(const Admixture& pca,
                              CountProfile<Alphabet>* p) const;
  // Rescales context specific profile with repeat penalty heuristic
  virtual void RescaleProfile(const matrix<double>& pp,
                              Profile<Alphabet>* profile) const;
  // Calculates posterior probability matrix for an input sequence.
  virtual void CalculatePosteriorsForSequence(const Sequence<Alphabet>& seq,
                                              matrix<double>* m) const;
  // Calculates posterior probability matrix for an input profile.
  virtual void CalculatePosteriorsForProfile(const CountProfile<Alphabet>& profile,
                                             matrix<double>* m) const;

 protected:
  // Needed to access names in templatized base class
  using ContextSpecificPseudocounts<Alphabet>::penalizer_;

 private:
  // Profile library with context profiles.
  const Hmm<Alphabet>& hmm_;
  // Needed to compute emission probabilities of context profiles.
  const MultEmission<Alphabet> emission_;

  DISALLOW_COPY_AND_ASSIGN(HMMPseudocounts);
};  // HMMPseudocounts

}  // namespace cs

#endif  // CS_HMM_PSEUDOCOUNTS_H_
