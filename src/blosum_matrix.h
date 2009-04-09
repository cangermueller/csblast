// Copyright 2009, Andreas Biegert

#ifndef SRC_BLOSUM_MATRIX_H_
#define SRC_BLOSUM_MATRIX_H_

#include <vector>

#include "amino_acid.h"
#include "substitution_matrix-inl.h"

namespace cs {

// BLOSUM family of substitution matrices for  class for substitution matrix
// classes.
class BlosumMatrix : public SubstitutionMatrix<AminoAcid> {
 public:
  enum Type {
    BLOSUM45 = 0,
    BLOSUM62 = 1,
    BLOSUM80 = 2
  };

  BlosumMatrix(Type matrix = BLOSUM62);
  virtual ~BlosumMatrix() { }

 private:
  // Initializes the matrix from target frequencies in raw data array.
  void init(const float* blosum_xx);
};

// Converts a BLOSUM matrix string to a BLOSUM matrix type.
BlosumMatrix::Type blosum_matrix_type_from_string(const std::string& s);

}  // namespace cs

#endif  // SRC_BLOSUM_MATRIX_H_
