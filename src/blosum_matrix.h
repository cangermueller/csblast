#ifndef CS_BLOSUM_MATRIX_H
#define CS_BLOSUM_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// BLOSUM family of substitution matrices for  class for substitution matrix classes.

#include <vector>

#include "amino_acid.h"
#include "exception.h"
#include "substitution_matrix.h"

namespace cs
{

class BlosumMatrix : public SubstitutionMatrix<AminoAcid>
{
  public:
    enum Type {
        BLOSUM45 = 0,
        BLOSUM62 = 1,
        BLOSUM80 = 2
    };

    BlosumMatrix(Type matrix = BLOSUM62);
    virtual ~BlosumMatrix() {}

  private:
    // Initializes the matrix from target frequencies in raw data array.
    void init(const float* blosum_xx);
};



// Converts a BLOSUM matrix string to a BLOSUM matrix type.
BlosumMatrix::Type blosum_matrix_type_from_string(const std::string& s);

}  // cs

#endif
