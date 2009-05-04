// Copyright 2009, Andreas Biegert

#ifndef SRC_CO_EMISSION_H_
#define SRC_CO_EMISSION_H_

#include <valarray>

#include "globals.h"
#include "profile.h"
#include "substitution_matrix.h"

namespace cs {

// Function object for calculation of co-emission score.
template< class Alphabet>
class CoEmission {
 public:
  // Constructs a co-emission object with background probabilities taken from
  // provided substitution matrix.
  CoEmission(const SubstitutionMatrix<Alphabet>* sm);

  ~CoEmission() {}

  // Calculates the normalized co-emission score of profile q and p over ncols
  // columns, starting at index qi in q and index pi in p. Both profiles must be
  // in lin-space.
  float operator() (const Profile<Alphabet>& q,
                     const Profile<Alphabet>& p,
                     int qi,
                     int pi,
                     int ncols) const;

 private:
  // Substition matrix needed for background frequencies
  const SubstitutionMatrix<Alphabet>* subst_matrix_;

  DISALLOW_COPY_AND_ASSIGN(CoEmission);
};  // class CoEmission

}  // namespace cs

#endif  // SRC_CO_EMISSION_H_
