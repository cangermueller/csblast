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
    enum Type {
        MATCH1_MISMATCH1 = 0,
        MATCH1_MISMATCH2 = 1,
        MATCH1_MISMATCH3 = 2,
        MATCH1_MISMATCH4 = 3,
        MATCH3_MISMATCH3 = 4,
        MATCH4_MISMATCH5 = 5
    };

    NucleotideMatrix(Type matrix = MATCH1_MISMATCH2);
    ~NucleotideMatrix() {}

private:
    // Disallow copy and assign.
    NucleotideMatrix(const NucleotideMatrix& other);
    NucleotideMatrix operator =(const NucleotideMatrix& other);

    // Initializes all matrix data members.
    void init(float match, float mismatch, const float* m);
};

}//cs

#endif
