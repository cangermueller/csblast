#ifndef CS_BLOSUM_MATRIX_H
#define CS_BLOSUM_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// BLOSUM family of substitution matrices for  class for substitution matrix classes.

#include <vector>

#include "substitution_matrix.h"

namespace cs
{

class BlosumMatrix : public SubstitutionMatrix
{
  public:
    enum Type {
        BLOSUM45 = 0,
        BLOSUM62 = 1,
        BLOSUM80 = 2
    };

    BlosumMatrix(Type matrix = BLOSUM62);
    ~BlosumMatrix() {}

private:
    // Disallow copy and assign.
    BlosumMatrix(const BlosumMatrix& other);
    BlosumMatrix operator =(const BlosumMatrix& other);

    // Set target frequencies to probabilities given in serialized matrix array
    void set_target_frequencies(const float* blosum_xx);
};

}//cs

#endif
