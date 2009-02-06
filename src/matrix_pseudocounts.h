#ifndef CS_MATRIX_PSEUDOCOUNTS_H
#define CS_MATRIX_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// An abstract base class for pseudocount methods.

#include "log.h"
#include "util.h"
#include "matrix.h"
#include "profile.h"
#include "counts_profile.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "substitution_matrix.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class Sequence;
template<class Alphabet_T>
class Profile;
template<class Alphabet_T>
class CountsProfile;

template<class Alphabet_T>
class MatrixPseudocounts : public Pseudocounts<Alphabet_T>
{
  public:
    MatrixPseudocounts(const SubstitutionMatrix<Alphabet_T>& m,
                       shared_ptr<Admixture> pca);
    ~MatrixPseudocounts() {}

    // Adds substitution matrix pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 Profile<Alphabet_T>& p) const;
    // Adds substitution matrix pseudocounts to alignment derived profile.
    virtual void add_to_profile(CountsProfile<Alphabet_T>& p) const;

  protected:
    // Needed to access names in templatized Profile base class
    using Pseudocounts<Alphabet_T>::admixture;

  private:
    // Disallow copy and assign
    MatrixPseudocounts(const MatrixPseudocounts&);
    void operator=(const MatrixPseudocounts&);

    // Substitution matrix with conditional probabilities for pseudocounts.
    const SubstitutionMatrix<Alphabet_T>& m_;
};  // MatrixPseudocounts



template<class Alphabet_T>
MatrixPseudocounts<Alphabet_T>::MatrixPseudocounts(const SubstitutionMatrix<Alphabet_T>& m,
                                                   shared_ptr<Admixture> pca)
        : Pseudocounts<Alphabet_T>(pca),
          m_(m)
{}

template<class Alphabet_T>
void MatrixPseudocounts<Alphabet_T>::add_to_sequence(const Sequence<Alphabet_T>& seq,
                                                     Profile<Alphabet_T>& p) const
{
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to sequence ...";
    LOG(DEBUG1) << p;
    if (seq.length() != p.num_cols())
        throw Exception("Cannot add substitution matrix pseudocounts: sequence and profile have different length!");

    // add substitution matrix pseudocounts
    float tau = admixture(1.0f);
    for(int i = 0; i < p.num_cols(); ++i) {
        for(int a = 0; a < p.alphabet_size(); ++a) {
            float pa = (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * m_.r(a, seq[i]);
            p[i][a] = p.logspace() ? log2(pa) : pa;
        }
    }
    normalize(p);
    LOG(DEBUG1) << p;
}

template<class Alphabet_T>
void MatrixPseudocounts<Alphabet_T>::add_to_profile(CountsProfile<Alphabet_T>& p) const
{
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to profile ...";
    LOG(DEBUG1) << p;
    const bool logspace = p.logspace();
    if (logspace) p.transform_to_linspace();

    // copy original frequencies to matrix f
    matrix<float> f(p.num_cols(), p.alphabet_size(), 0.0f);
    for (int i = 0; i < p.num_cols(); ++i)
        for(int a = 0; a < p.alphabet_size(); ++a)
            f[i][a] = p[i][a];

    // add substitution matrix pseudocounts
    for(int i = 0; i < p.num_cols(); ++i) {
        float tau = admixture(p.neff(i));
        for(int a = 0; a < p.alphabet_size(); ++a) {
            float sum = 0.0f;
            for(int b = 0; b < p.alphabet_size(); ++b)
                sum += m_.r(a,b) * f[i][b];
            p[i][a] = (1.0f - tau) * f[i][a] + tau * sum;
        }
    }

    if (logspace) p.transform_to_logspace();
    normalize(p);
    LOG(DEBUG1) << p;
}

}  // cs

#endif
