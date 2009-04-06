#ifndef CS_MATRIX_PSEUDOCOUNTS_H
#define CS_MATRIX_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of substitution matrix pseudocounts.

#include "counts_profile.h"
#include "log.h"
#include "matrix.h"
#include "profile.h"
#include "pseudocounts.h"
#include "sequence.h"
#include "substitution_matrix.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class MatrixPseudocounts : public Pseudocounts<Alphabet_T>
{
  public:
    MatrixPseudocounts(const SubstitutionMatrix<Alphabet_T>* m);
    ~MatrixPseudocounts() {}

    // Adds substitution matrix pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 const Admixture& pca,
                                 Profile<Alphabet_T>* profile) const;
    // Adds substitution matrix pseudocounts to alignment derived profile.
    virtual void add_to_profile(const Admixture& pca, CountsProfile<Alphabet_T>* profile) const;

  private:
    // Disallow copy and assign
    MatrixPseudocounts(const MatrixPseudocounts&);
    void operator=(const MatrixPseudocounts&);

    // Substitution matrix with conditional probabilities for pseudocounts.
    const SubstitutionMatrix<Alphabet_T>* m_;
};  // MatrixPseudocounts



template<class Alphabet_T>
inline MatrixPseudocounts<Alphabet_T>::MatrixPseudocounts(const SubstitutionMatrix<Alphabet_T>* m)
        : m_(m)
{ }

template<class Alphabet_T>
void MatrixPseudocounts<Alphabet_T>::add_to_sequence(const Sequence<Alphabet_T>& seq,
                                                     const Admixture& pca,
                                                     Profile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding substitution matrix pseudocounts to sequence ...";
    if (seq.length() != profile->num_cols())
        throw Exception("Cannot add substitution matrix pseudocounts: sequence and profile have different length!");

    // add substitution matrix pseudocounts
    Profile<Alphabet_T>& p = *profile;
    float tau = pca(1.0f);
    for(int i = 0; i < p.num_cols(); ++i) {
        for(int a = 0; a < p.alphabet_size(); ++a) {
            float pa = (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * m_->r(a, seq[i]);
            p[i][a] = p.logspace() ? log2(pa) : pa;
        }
    }
    normalize(profile);
    LOG(DEBUG2) << *profile;
}

template<class Alphabet_T>
void MatrixPseudocounts<Alphabet_T>::add_to_profile(const Admixture& pca,
                                                    CountsProfile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding substitution matrix pseudocounts to profile ...";

    CountsProfile<Alphabet_T>& p = *profile;
    const bool logspace = p.logspace();
    if (logspace) p.transform_to_linspace();

    // copy original frequencies to matrix f
    matrix<float> f(p.num_cols(), p.alphabet_size(), 0.0f);
    for (int i = 0; i < p.num_cols(); ++i)
        for(int a = 0; a < p.alphabet_size(); ++a)
            f[i][a] = p[i][a];

    // add substitution matrix pseudocounts
    for(int i = 0; i < p.num_cols(); ++i) {
        float tau = pca(p.neff(i));
        for(int a = 0; a < p.alphabet_size(); ++a) {
            float sum = 0.0f;
            for(int b = 0; b < p.alphabet_size(); ++b)
                sum += m_->r(a,b) * f[i][b];
            p[i][a] = (1.0f - tau) * f[i][a] + tau * sum;
        }
    }

    if (logspace) p.transform_to_logspace();
    normalize(profile);
    LOG(DEBUG2) << *profile;
}

}  // cs

#endif
