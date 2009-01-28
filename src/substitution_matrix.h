#ifndef CS_SUBSTITUTION_MATRIX_H
#define CS_SUBSTITUTION_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for substitution matrix classes.

#include <cmath>

#include <iostream>
#include <iomanip>
#include <vector>

#include "matrix.h"
#include "sequence_alphabet.h"

namespace cs
{

class SubstitutionMatrix
{
  public:
        // To be used by derived classes.
    SubstitutionMatrix(const SequenceAlphabet* alphabet);
    ~SubstitutionMatrix() {}

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

    friend std::ostream& operator<< (std::ostream& out, const SubstitutionMatrix& m);

  protected:
    // Initializes all matrix data members.
    virtual void init() = 0;
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
    // Sequence alphabet over which the matrix records substitution scores.
    const SequenceAlphabet* alphabet_;
};

// Prints the substitution matrix in human readable format to stream.
inline std::ostream& operator<< (std::ostream& out, const SubstitutionMatrix& m)
{
    out << "Background frequencies:\n";
    out << *m.alphabet_ << std::endl;
    for (int a = 0; a < m.size_; ++a)
        out << std::fixed << std::setprecision(1) << 100 * m.f_[a] << "\t";
    out << "\nSubstitution matrix log2( P(a,b) / p(a)*p(b) ) (in bits):\n";
    out << *m.alphabet_ << std::endl;
    for (int b = 0; b < m.size_; ++b) {
        for (int a = 0; a < m.size_; ++a)
            out << std::fixed << std::setprecision(1) << std::showpos << m.s_[a][b] << std::noshowpos << "\t";
        out << std::endl;
    }
    out << "Probability matrix P(a,b) (in %):\n";
    out << *m.alphabet_ << std::endl;
    for (int b = 0; b < m.size_; ++b) {
        for (int a = 0; a < m.size_; ++a)
            out << std::fixed << std::setprecision(1) << 100 * m.p_[b][a] << "\t";
        out << std::endl;
    }
    out << "Matrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
    out << *m.alphabet_ << std::endl;
    for (int b = 0; b < m.size_; ++b) {
        for (int a = 0; a < m.size_; ++a)
            out << std::fixed << std::setprecision(1) << 100 * m.r_[b][a] << "\t";
        out << std::endl;
    }
    return out;
}

}  // cs

#endif
