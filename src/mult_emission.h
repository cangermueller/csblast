// Copyright 2009, Andreas Biegert

#ifndef SRC_MULT_EMISSION_H_
#define SRC_MULT_EMISSION_H_

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
                     const CountProfile<Alphabet>& counts_profile,
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

#endif  // SRC_MULT_EMISSION_H_
