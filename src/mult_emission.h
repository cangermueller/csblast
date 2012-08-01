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

#ifndef CS_MULT_EMISSION_H_
#define CS_MULT_EMISSION_H_

#include <valarray>

#include "globals.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "sequence-inl.h"

namespace cs {

// Function object for calculation of multinomial emission probabilities for
// context profiles.
template< class Alphabet>
class MultEmission {
 public:
  // Constructs an emission object with positional window weights.
  MultEmission(int num_cols, float weight_center = 1.0f, float weight_decay = 1.0f);

  ~MultEmission() {}

  // Calculates the log emission probability of profile window centered at given
  // index.
  double operator() (const ContextProfile<Alphabet>& profile,
                     const CountProfile<Alphabet>& count_profile,
                     int index) const;
  // Calculates the log emission probability of sequence window centered at given
  // index.
  double operator() (const ContextProfile<Alphabet>& profile,
                     const Sequence<Alphabet>& seq,
                     int index) const;
  // Calculates the sum of positional weights.
  float SumWeights() const { return weights_.sum(); }
  // Initializes positional window weights.
  void InitWeights(float weight_center, float weight_decay);

 private:
  // Number of columns in context profiles.
  int num_cols_;
  // Index of central column in context profiles.
  int center_;
  // Positional window weights
  std::valarray<float> weights_;

  DISALLOW_COPY_AND_ASSIGN(MultEmission);
};  // class MultEmission

}  // namespace cs

#endif  // CS_MULT_EMISSION_H_
