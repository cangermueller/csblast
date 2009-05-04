// Copyright 2009, Andreas Biegert

#ifndef SRC_CO_EMISSION_INL_H_
#define SRC_CO_EMISSION_INL_H_

#include "co_emission.h"

#include <cassert>

#include "profile-inl.h"
#include "substitution_matrix.h"
#include "log.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
CoEmission<Alphabet>::CoEmission(const SubstitutionMatrix<Alphabet>* sm)
    : subst_matrix_(sm) {}

template<class Alphabet>
inline float CoEmission<Alphabet>::operator() (const Profile<Alphabet>& q,
                                               const Profile<Alphabet>& p,
                                               int qi,
                                               int pi,
                                               int ncols) const {
  assert(!q.logspace());
  assert(!p.logspace());

  const int alphabet_size = Alphabet::instance().size();
  float rv = 0.0f;
  int i = qi;
  int j = pi;

  for (int k = 0; k < ncols; ++k) {
    float sum = 0.0f;
    for (int a = 0; a < alphabet_size; ++a)
      sum += q[i][a] * p[j][a] / subst_matrix_->f(a);

    rv += fast_log2(sum);
    ++i;
    ++j;
  }

  rv /= ncols;
  return rv;
}

}  // namespace cs

#endif  // SRC_CO_EMISSION_INL_H_
