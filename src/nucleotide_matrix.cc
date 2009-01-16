/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleotide_matrix.h"

#include <iostream>

#include "nucleotide_alphabet.h"
#include "matrix.h"
#include "substitution_matrix.h"
#include "util.h"

namespace
{

const float g_match1_mismatch2[] = {
    //A    C      G      T
    0.51948052,
    -0.05194805, 0.51948052,
    -0.05194805, -0.05194805, 0.51948052,
    -0.05194805, -0.05194805, -0.05194805,  0.51948052 };

}  // namnespace

namespace cs
{

NucleotideMatrix::NucleotideMatrix(Type matrix)
        : SubstitutionMatrix(NucleotideAlphabet::instance()->size() - 1)
{
    switch (matrix) {
        case MATCH1_MISMATCH2:
            init(1, -2, g_match1_mismatch2);
            break;
        default:
            throw MyException("Unsupported BLOSUM matrix!");
    }
}

void NucleotideMatrix::init(float match, float mismatch, const float* m)
{
    // Read raw data vector
    Matrix<float> m_inv(size_, size_, 0.0f);
    int n = 0;
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b <= a; ++b, ++n)
            m_inv(a,b) = m[n];
    // Add uppper right matrix part
    for (int a = 0; a < size_-1; ++a)
        for (int b = a+1; b < size_; ++b)
            m_inv(a,b) = m_inv(b,a);

    // Calculate background frequencies
    for (int a = 0; a < size_; ++a) {
        f_[a] = 0.0f;
        for (int b = 0; b < size_; ++b) f_[a] += m_inv(a,b);
    }
    normalize_to_one(&f_[0], size_);

    // Fille substitution matrix
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            s_(a,b) = a==b ? match : mismatch;

    init_from_substitution_matrix_and_background_frequencies();
}

}  // cs
