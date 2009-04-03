#ifndef CS_LIBRARY_PSEUDOCOUNTS_H
#define CS_LIBRARY_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of context-specific pseudocounts calculated from a profile
// library.

#include "counts_profile.h"
#include "emitter.h"
#include "log.h"
#include "matrix.h"
#include "profile.h"
#include "profile_library.h"
#include "pseudocounts.h"
#include "sequence.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class LibraryPseudocounts : public Pseudocounts<Alphabet_T>
{
  public:
    LibraryPseudocounts(const ProfileLibrary<Alphabet_T>* lib);
    ~LibraryPseudocounts() { }

    // Adds context-specific pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 const Admixture& pca,
                                 Profile<Alphabet_T>* profile) const;
    // Adds context-specific pseudocounts to alignment derived profile.
    virtual void add_to_profile(const Admixture& pca, CountsProfile<Alphabet_T>* p) const;
    // Sets substiution matrix to be used for conditional probabilities.
    void set_matrix(const ProfileLibrary<Alphabet_T>* lib) { lib_ = lib; }

  private:
    // Disallow copy and assign
    LibraryPseudocounts(const LibraryPseudocounts&);
    void operator=(const LibraryPseudocounts&);

    // Profile library with context profiles.
    const ProfileLibrary<Alphabet_T>* lib_;
};  // LibraryPseudocounts



template<class Alphabet_T>
LibraryPseudocounts<Alphabet_T>::LibraryPseudocounts(const ProfileLibrary<Alphabet_T>* lib)
        : lib_(lib)
{ }

template<class Alphabet_T>
void LibraryPseudocounts<Alphabet_T>::add_to_sequence(const Sequence<Alphabet_T>& seq,
                                                      const Admixture& pca,
                                                      Profile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding context-specific, library derived pseudocounts to sequence ...";
    LOG(DEBUG2) << *profile;
    if (seq.length() != profile->num_cols())
        throw Exception("Cannot add context-specific pseudocounts: sequence and profile have different length!");

    float tau = pca(1.0f);
    CountsProfile<Alphabet_T>& p = *profile;
    for(int i = 0; i < p.num_cols(); ++i) {
        for(int a = 0; a < p.alphabet_size(); ++a) {
            float pa = (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * m_->r(a, seq[i]);
            p[i][a] = p.logspace() ? log2(pa) : pa;
        }
    }
    normalize(p);
    LOG(DEBUG2) << p;
}

template<class Alphabet_T>
void LibraryPseudocounts<Alphabet_T>::add_to_profile(const Admixture& pca, CountsProfile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding substitution matrix pseudocounts to profile ...";
    LOG(DEBUG2) << p;
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
    normalize(p);
    LOG(DEBUG2) << p;
}

}  // cs

#endif
