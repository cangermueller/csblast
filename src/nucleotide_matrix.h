#ifndef CS_NUCLEOTIDE_MATRIX_H
#define CS_NUCLEOTIDE_MATRIX_H
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

class NucleotideMatrix : public SubstitutionMatrix
{
  public:
    NucleotideMatrix(float match, float mismatch);
    ~NucleotideMatrix() {}

private:
    // Disallow copy and assign.
    NucleotideMatrix(const NucleotideMatrix& other);
    NucleotideMatrix operator =(const NucleotideMatrix& other);

    // Initializes all matrix data members.
    void init(float match, float mismatch);
};

}//cs

#endif
