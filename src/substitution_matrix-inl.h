// Copyright 2009, Andreas Biegert

#ifndef SRC_SUBSTITUTION_MATRIX_INL_H_
#define SRC_SUBSTITUTION_MATRIX_INL_H_

#include "substitution_matrix.h"

#include <iomanip>
#include <iostream>

#include "log.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
SubstitutionMatrix<Alphabet>::SubstitutionMatrix()
    : size_(Alphabet::instance().size()),
      p_(size_, size_, 0.0f),
      s_(size_, size_, 0.0f),
      r_(size_, size_, 0.0f),
      f_(size_, 0.0f) {}

template<class Alphabet>
SubstitutionMatrix<Alphabet>::~SubstitutionMatrix() {}

template<class Alphabet>
void SubstitutionMatrix<Alphabet>::init_from_target_frequencies() {
  // Check transition probability matrix, renormalize P
  float sumab = 0.0f;
  for (int a = 0; a < size_; a++)
    for (int b = 0; b < size_; ++b) sumab += p_[a][b];
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b) p_[a][b] /= sumab;

  // Calculate background frequencies
  for (int a = 0; a < size_; ++a) {
    f_[a] = 0.0f;
    for (int b = 0; b < size_; ++b) f_[a] += p_[a][b];
  }
  normalize_to_one(&f_[0], size_);

  // Precompute matrix R for amino acid pseudocounts:
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b)
      r_[a][b] = p_[a][b] / f_[b]; // R[a][b] = P(a|b)

  // Calculate scoring matrix
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b)
      // S[a][b] = log2(P(a,b) / P(a)*P(b))
      s_[a][b] = log2(r_[a][b] / f_[a]);

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void SubstitutionMatrix<Alphabet>::init_from_scores_and_background_frequencies() {
  // Calculate target frequencies
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b)
      p_[a][b] = fast_pow2(s_[a][b]) * f_[a] * f_[b];

  // Check transition probability matrix, renormalize P
  float sumab = 0.0f;
  for (int a = 0; a < size_; a++)
    for (int b = 0; b < size_; ++b) sumab += p_[a][b];
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b) p_[a][b] /= sumab;

  // Precompute matrix R for amino acid pseudocounts:
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b)
      r_[a][b] = p_[a][b] / f_[b]; // R[a][b] = P(a|b)

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void SubstitutionMatrix<Alphabet>::print(std::ostream& out) const {
  std::ios_base::fmtflags flags = out.flags();  // save flags

  out << "Background frequencies:\n";
  out << Alphabet::instance() << std::endl;
  for (int a = 0; a < size_; ++a)
    out << std::fixed << std::setprecision(1) << 100 * f_[a] << "\t";
  out << "\nSubstitution trix log2( P(a,b) / p(a)*p(b) ) (in bits):\n";
  out << Alphabet::instance() << std::endl;
  for (int b = 0; b < size_; ++b) {
    for (int a = 0; a < size_; ++a)
      out << std::fixed << std::setprecision(1) << std::showpos
          << s_[a][b] << std::noshowpos << "\t";
    out << std::endl;
  }
  out << "Probability matrix P(a,b) (in %):\n";
  out << Alphabet::instance() << std::endl;
  for (int b = 0; b < size_; ++b) {
    for (int a = 0; a < size_; ++a)
      out << std::fixed << std::setprecision(1) << 100 * p_[b][a] << "\t";
    out << std::endl;
  }
  out << "Matrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
  out << Alphabet::instance() << std::endl;
  for (int b = 0; b < size_; ++b) {
    for (int a = 0; a < size_; ++a)
      out << std::fixed << std::setprecision(1) << 100 * r_[b][a] << "\t";
    out << std::endl;
  }

  out.flags(flags);
}

}  // namespace cs

#endif  // SRC_SUBSTITUTION_MATRIX_INL_H_
