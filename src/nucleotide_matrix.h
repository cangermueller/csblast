// Copyright 2009, Andreas Biegert

#ifndef SRC_NUCLEOTIDE_MATRIX_H_
#define SRC_NUCLEOTIDE_MATRIX_H_

#include "nucleotide.h"
#include "substitution_matrix-inl.h"

namespace cs {

// BLOSUM family of substitution matrices for  class for substitution matrix
// classes.
class NucleotideMatrix : public SubstitutionMatrix<Nucleotide> {
 public:
  NucleotideMatrix(float match, float mismatch) {
    Init(match, mismatch);
  }
  virtual ~NucleotideMatrix() {}

 private:
  // Initializes all matrix data members.
  virtual void Init(float match, float mismatch);
};

}  // namespace cs

#endif  // SRC_NUCLEOTIDE_MATRIX_H_
