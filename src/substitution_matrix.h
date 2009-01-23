#ifndef CS_SUBSTITUTION_MATRIX_H
#define CS_SUBSTITUTION_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for substitution matrix classes.

#include <cmath>

#include <vector>

#include "matrix.h"

namespace cs
{

class SubstitutionMatrix
{
  public:
    // Return substitution score S(a,b).
    float s(int a, int b) const { return s_[a][b]; }
    // Returns joint probability P(a,b).
    float p(int a, int b) const { return p_[a][b]; }
    // Returns probability P(a|b) "a given b".
    float r(int a, int b) const { return r_[a][b]; }
    // Returns background frequency of a.
    float f(int a) const { return f_[a]; }

    // Returns the size of the substitution matrix.
    int size() const { return size_; }

protected:
    // To be used by derived classes.
    SubstitutionMatrix(int size);
    ~SubstitutionMatrix() {}

    // Initializes other matrix data members from matrix P.
    void init_from_target_frequencies();
    // Initializes other matrix data members from substitution matrix S and background frequencies.
    void init_from_substitution_matrix_and_background_frequencies();

    // Size of the matrix.
    const int size_;
    // Target frequency matrix P(a,b).
    matrix<float> p_;
    // Substitution matrix S(a,b).
    matrix<float> s_;
    // "a given b" probability matrix.
    matrix<float> r_;
    // Background frequencies of alphabet.
    std::vector<float> f_;

private:
    // Disallow copy and assign.
    SubstitutionMatrix(const SubstitutionMatrix& other);
    SubstitutionMatrix operator =(const SubstitutionMatrix& other);

    void print_debug() const;
};

}//cs

#endif
