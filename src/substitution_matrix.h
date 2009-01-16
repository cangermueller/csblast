#ifndef CS_SUBSTITUTION_MATRIX_H
#define CS_SUBSTITUTION_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for substitution matrix classes.

#include <vector>

#include "matrix.h"

namespace cs
{

class SubstitutionMatrix
{
  public:
    // Return substitution score s(a,b).
    float s(int a, int b) const { return s_(a,b); }
    // Returns joint probability p(a,b).
    float p(int a, int b) const { return p_(a,b); }
    // Returns probability p(a|b) "a given b".
    float r(int a, int b) const { return r_(a,b); }
    // Returns background frequency of a.
    float f(int a) const { return f_[a]; }

    // Returns the size of the substitution matrix.
    int size() const { return size_; }

protected:
    // Size of the matrix.
    int size_;
    // Substitution matrix s(a,b).
    Matrix<float> s_;
    // Joint probability matrix p(a,b).
    Matrix<float> p_;
    // "a given b" probability matrix.
    Matrix<float> r_;
    // Background frequencies of alphabet.
    std::vector<float> f_;

private:
    // Disallow copy and assign.
    SubstitutionMatrix(const SubstitutionMatrix& other);
    SubstitutionMatrix operator =(const SubstitutionMatrix& other);
};

}//cs

#endif
