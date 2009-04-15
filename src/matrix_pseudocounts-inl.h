// Copyright 2009, Andreas Biegert

#ifndef SRC_MATRIX_PSEUDOCOUNTS_INL_H_
#define SRC_MATRIX_PSEUDOCOUNTS_INL_H_

#include "matrix_pseudocounts.h"

#include "count_profile-inl.h"
#include "log.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
inline MatrixPseudocounts<Alphabet>::MatrixPseudocounts(
    const SubstitutionMatrix<Alphabet>* m)
    : m_(m) {}

template<class Alphabet>
void MatrixPseudocounts<Alphabet>::add_to_sequence(
    const Sequence<Alphabet>& seq,
    const Admixture& pca,
    Profile<Alphabet>* profile) const {
  LOG(DEBUG2) << "Adding substitution matrix pseudocounts to sequence ...";

  if (seq.length() != profile->num_cols())
    throw Exception("Cannot add substitution matrix pseudocounts: "
                    "sequence and profile have different length!");

  // add substitution matrix pseudocounts
  Profile<Alphabet>& p = *profile;
  float tau = pca(1.0f);
  for(int i = 0; i < p.num_cols(); ++i) {
    for(int a = 0; a < p.alphabet_size(); ++a) {
      float pa =
        (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f)
        + tau * m_->r(a, seq[i]);
      p[i][a] = p.logspace() ? fast_log2(pa) : pa;
    }
  }
  normalize(profile);

  LOG(DEBUG2) << *profile;
}

template<class Alphabet>
void MatrixPseudocounts<Alphabet>::add_to_profile(
    const Admixture& pca,
    CountProfile<Alphabet>* profile) const {
  LOG(DEBUG2) << "Adding substitution matrix pseudocounts to profile ...";

  CountProfile<Alphabet>& p = *profile;
  const bool logspace = p.logspace();
  if (logspace) p.transform_to_linspace();

  // copy original frequencies to matrix f
  matrix<float> f(p.num_cols(), p.alphabet_size(), 0.0f);
  for (int i = 0; i < p.num_cols(); ++i)
    for(int a = 0; a < p.alphabet_size(); ++a)
      f[i][a] = p[i][a];

  // add substitution matrix pseudocounts
  for(int i = 0; i < p.num_cols(); ++i) {
    float tau = pca(p.neff(i));
    for(int a = 0; a < p.alphabet_size(); ++a) {
      float sum = 0.0f;
      for(int b = 0; b < p.alphabet_size(); ++b)
        sum += m_->r(a,b) * f[i][b];
      p[i][a] = (1.0f - tau) * f[i][a] + tau * sum;
    }
  }

  if (logspace) p.transform_to_logspace();
  normalize(profile);

  LOG(DEBUG2) << *profile;
}

}  // namespace cs

#endif  // SRC_MATRIX_PSEUDOCOUNTS_INL_H_
