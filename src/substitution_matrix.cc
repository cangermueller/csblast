/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "substitution_matrix.h"

#include <cstdio>

#include <iostream>

#include "util.h"

namespace
{

const bool kDebug = false;

} // namespace

namespace cs
{

SubstitutionMatrix::SubstitutionMatrix(int size)
  : size_(size),
    p_(size, size, 0.0f),
    s_(size, size, 0.0f),
    r_(size, size, 0.0f),
    f_(size, 0.0f)
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

    if (kDebug) print_debug();
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

    if (kDebug) print_debug();
}

void SubstitutionMatrix::print_debug() const
{
    std::cerr << std::endl << "Background frequencies:\nf[] ";
    for (int a = 0; a < size_; ++a) fprintf(stderr, "%4.1f ", 100 * f_[a]);
    std::cerr << std::endl << "\nSubstitution matrix log2( P(a,b) / p(a)*p(b) ) (in bits):\n";
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)  fprintf(stderr, "%4.1f ", s_[a][b]);
        std::cerr << std::endl;
    }
    std::cerr << std::endl << "Probability matrix P(a,b) (in %):\n";
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)  fprintf(stderr, "%4.1f ", 100 * p_[b][a]);
        std::cerr << std::endl;
    }
    std::cerr << std::endl << "Matrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)  fprintf(stderr, "%4.1f ", 100 * r_[b][a]);
        std::cerr << std::endl;
    }
}

}  // cs
