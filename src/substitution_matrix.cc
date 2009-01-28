/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "substitution_matrix.h"

#include <cstdio>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "log.h"
#include "util.h"
#include "sequence_alphabet.h"

namespace cs
{

SubstitutionMatrix::SubstitutionMatrix(const SequenceAlphabet* alphabet)
        : size_(alphabet->size()),
          p_(size_, size_, 0.0f),
          s_(size_, size_, 0.0f),
          r_(size_, size_, 0.0f),
          f_(size_, 0.0f),
          alphabet_(alphabet)
{}

void SubstitutionMatrix::init_from_target_frequencies()
{
    // Check transition probability matrix, renormalize P
    float sumab = 0.0f;
    for (int a = 0; a < size_; a++)
        for (int b = 0; b < size_; ++b) sumab += p_[a][b];
    for (int a = 0; a < size_; ++a)
         for (int b = 0; b < size_; ++b) p_[a][b] /= sumab;

    // Calculate background frequencies
    for (int a = 0; a < size_; ++a) {
        f_[a] = 0.0f;
        for (int b = 0; b < size_; ++b) f_[a] += p_[a][b];
    }
    normalize_to_one(&f_[0], size_);

    // Precompute matrix R for amino acid pseudocounts:
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            r_[a][b] = p_[a][b] / f_[b]; // R[a][b] = P(a|b)

    // Calculate scoring matrix
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            s_[a][b] = log2(r_[a][b] / f_[a]); // S[a][b] = log2(P(a,b) / P(a)*P(b))

    LOG(DEBUG1) << *this;
}

void SubstitutionMatrix::init_from_substitution_matrix_and_background_frequencies()
{
    // Calculate target frequencies
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            p_[a][b] = pow(2.0, s_[a][b]) * f_[a] * f_[b];

    // Check transition probability matrix, renormalize P
    float sumab = 0.0f;
    for (int a = 0; a < size_; a++)
        for (int b = 0; b < size_; ++b) sumab += p_[a][b];
    for (int a = 0; a < size_; ++a)
         for (int b = 0; b < size_; ++b) p_[a][b] /= sumab;

    // Precompute matrix R for amino acid pseudocounts:
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            r_[a][b] = p_[a][b] / f_[b]; // R[a][b] = P(a|b)

    LOG(DEBUG1) << *this;
}

}  // cs
