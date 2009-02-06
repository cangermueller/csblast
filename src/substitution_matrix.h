#ifndef CS_SUBSTITUTION_MATRIX_H
#define CS_SUBSTITUTION_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for substitution matrix classes.

#include <cmath>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "log.h"
#include "matrix.h"
#include "util.h"

namespace cs
{

template<class Alphabet_T>
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

    // Prints the substitution matrix in human readable format to stream.
    friend std::ostream& operator<< (std::ostream& out, const SubstitutionMatrix& m)
    {
        m.print(out);
        return out;
    }

  protected:
    // To be used by derived classes.
    SubstitutionMatrix();
    virtual ~SubstitutionMatrix() = 0;

    // Initializes the other matrix data members from target frequencies.
    void init_from_target_frequencies();
    // Initializes the other matrix data members from score matrix and background frequencies.
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
    // Prints the substitution matrix in human-readable format to output stream.
    virtual void print(std::ostream& out) const;
};



template<class Alphabet_T>
SubstitutionMatrix<Alphabet_T>::SubstitutionMatrix()
        : size_(Alphabet_T::instance().size()),
          p_(size_, size_, 0.0f),
          s_(size_, size_, 0.0f),
          r_(size_, size_, 0.0f),
          f_(size_, 0.0f)
{}

template<class Alphabet_T>
SubstitutionMatrix<Alphabet_T>::~SubstitutionMatrix()
{}

template<class Alphabet_T>
void SubstitutionMatrix<Alphabet_T>::init_from_target_frequencies()
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

template<class Alphabet_T>
void SubstitutionMatrix<Alphabet_T>::init_from_substitution_matrix_and_background_frequencies()
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

template<class Alphabet_T>
void SubstitutionMatrix<Alphabet_T>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save flags

    out << "Background frequencies:\n";
    out << Alphabet_T::instance() << std::endl;
    for (int a = 0; a < size_; ++a)
        out << std::fixed << std::setprecision(1) << 100 * f_[a] << "\t";
    out << "\nSubstitution trix log2( P(a,b) / p(a)*p(b) ) (in bits):\n";
    out << Alphabet_T::instance() << std::endl;
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)
            out << std::fixed << std::setprecision(1) << std::showpos << s_[a][b] << std::noshowpos << "\t";
        out << std::endl;
    }
    out << "Probability matrix P(a,b) (in %):\n";
    out << Alphabet_T::instance() << std::endl;
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)
            out << std::fixed << std::setprecision(1) << 100 * p_[b][a] << "\t";
        out << std::endl;
    }
    out << "Matrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
    out << Alphabet_T::instance() << std::endl;
    for (int b = 0; b < size_; ++b) {
        for (int a = 0; a < size_; ++a)
            out << std::fixed << std::setprecision(1) << 100 * r_[b][a] << "\t";
        out << std::endl;
    }

    out.flags(flags);
}

}  // cs

#endif
