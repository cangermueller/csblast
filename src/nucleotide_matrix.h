#ifndef CS_NUCLEOTIDE_MATRIX_H
#define CS_NUCLEOTIDE_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// BLOSUM family of substitution matrices for  class for substitution matrix classes.

#include "nucleotide.h"
#include "substitution_matrix.h"

namespace cs
{

class NucleotideMatrix : public SubstitutionMatrix<Nucleotide>
{
  public:
    NucleotideMatrix(float match, float mismatch);
    virtual ~NucleotideMatrix() {}

  private:
    // Initializes all matrix data members.
    virtual void init(float match, float mismatch);
};

}  // cs

#endif
