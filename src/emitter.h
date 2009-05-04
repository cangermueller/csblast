// Copyright 2009, Andreas Biegert

#ifndef SRC_EMITTER_H_
#define SRC_EMITTER_H_

#include <valarray>

#include "globals.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation for computation of emitter probabilities for profiles.
template< class Alphabet>
class Emitter {
 public:
  // Constructs an emitter with positional window weights.
  Emitter(int num_cols, float weight_center = 1.0f, float weight_decay = 1.0f);

  ~Emitter() {}

  // Calculates the log emission probability of profile window centered at given
  // index.
  double operator() (const ContextProfile<Alphabet>& profile,
                     const CountProfile<Alphabet>& counts,
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

  DISALLOW_COPY_AND_ASSIGN(Emitter);
};  // class Emitter

}  // namespace cs

#endif  // SRC_EMITTER_H_
