// Copyright 2009, Andreas Biegert

#ifndef SRC_MATRIX_PSEUDOCOUNTS_H_
#define SRC_MATRIX_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "globals.h"
#include "log.h"
#include "matrix.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Substitution matrix pseudocounts factory.
template<class Alphabet>
class MatrixPseudocounts : public Pseudocounts<Alphabet> {
 public:
  MatrixPseudocounts(const SubstitutionMatrix<Alphabet>* m);
  ~MatrixPseudocounts() {}

  // Adds substitution matrix pseudocounts to sequence and stores resulting
  // frequencies in given profile.
  virtual void add_to_sequence(const Sequence<Alphabet>& seq,
                               const Admixture& pca,
                               Profile<Alphabet>* profile) const;
  // Adds substitution matrix pseudocounts to alignment derived profile.
  virtual void add_to_profile(const Admixture& pca,
                              CountProfile<Alphabet>* profile) const;

 private:
  // Substitution matrix with conditional probabilities for pseudocounts.
  const SubstitutionMatrix<Alphabet>* m_;

  DISALLOW_COPY_AND_ASSIGN(MatrixPseudocounts);
};  // MatrixPseudocounts

}  // namespace cs

#endif  // SRC_MATRIX_PSEUDOCOUNTS_H_
