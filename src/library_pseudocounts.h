#ifndef CS_LIBRARY_PSEUDOCOUNTS_H
#define CS_LIBRARY_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of context-specific pseudocounts calculated from a library
// of context profiles.

#include <cassert>
#include <cmath>

#include <valarray>

#include "count_profile-inl.h"
#include "emitter.h"
#include "log.h"
#include "matrix.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class LibraryPseudocounts : public Pseudocounts<Alphabet_T>
{
  public:
    LibraryPseudocounts(const ProfileLibrary<Alphabet_T>* lib, const EmissionParams& params);
    ~LibraryPseudocounts() { }

    // Adds context-specific pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 const Admixture& pca,
                                 Profile<Alphabet_T>* profile) const;
    // Adds context-specific pseudocounts to alignment derived profile.
    virtual void add_to_profile(const Admixture& pca, CountProfile<Alphabet_T>* p) const;

  private:
    // Disallow copy and assign
    LibraryPseudocounts(const LibraryPseudocounts&);
    void operator=(const LibraryPseudocounts&);

    // Profile library with context profiles.
    const ProfileLibrary<Alphabet_T>& lib_;
    // Needed to compute emission probabilities of context profiles.
    const Emitter<Alphabet_T> emitter_;
};  // LibraryPseudocounts



template<class Alphabet_T>
inline LibraryPseudocounts<Alphabet_T>::LibraryPseudocounts(const ProfileLibrary<Alphabet_T>* lib,
                                                            const EmissionParams& params)
        : lib_(*lib),
          emitter_(lib->num_cols(), params)
{
    assert(lib_.logspace());
}

template<class Alphabet_T>
void LibraryPseudocounts<Alphabet_T>::add_to_sequence(const Sequence<Alphabet_T>& seq,
                                                      const Admixture& pca,
                                                      Profile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding context-specific, library-derived pseudocounts to sequence ...";
    LOG(DEBUG2) << *profile;
    if (seq.length() != profile->num_cols())
        throw Exception("Cannot add context-specific pseudocounts: sequence and profile have different length!");

    const int length        = seq.length();
    const int num_profiles  = lib_.num_profiles();
    const int alphabet_size = Alphabet_T::instance().size();
    const int center        = lib_.center();
    const float tau         = pca(1.0f);  // number of effective seqs is one

    std::valarray<float> prob(0.0f, num_profiles);  // profile probabilities P(p_k|X_i) at position i
    std::valarray<float> pc(0.0f, alphabet_size);   // pseudocount vector P(a|X_i) at position i
    Profile<Alphabet_T>& p = *profile;              // output profile with pseudocounts

    for (int i = 0; i < length; ++i) {
        // calculate profile probabilities P(p_k|X_i)
        for (int k = 0; k < num_profiles; ++k) {
            prob[k] = lib_[k].prior() * pow(2.0, emitter_(lib_[k], seq, i));
        }
        prob /= prob.sum();  // normalization

        // calculate pseudocount vector P(a|X_i)
        for(int a = 0; a < alphabet_size; ++a) {
            pc[a] = 0.0f;
            for (int k = 0; k < num_profiles; ++k) {
                pc[a] += prob[k] * pow(2.0, lib_[k][center][a]);
            }
        }
        pc /= pc.sum();  // normalization

        // add pseudocounts to sequence by storing probabilities in output profile
        for(int a = 0; a < alphabet_size; ++a) {
            float pa = (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * pc[a];
            p[i][a] = p.logspace() ? log2(pa) : pa;
        }
    }

    normalize(profile);
    LOG(DEBUG2) << *profile;
}

template<class Alphabet_T>
void LibraryPseudocounts<Alphabet_T>::add_to_profile(const Admixture& pca, CountProfile<Alphabet_T>* profile) const
{
    assert(!profile->has_counts());
    assert(!profile->logspace());

    LOG(DEBUG2) << "Adding context-specific, library-derived pseudocounts to profile ...";
    LOG(DEBUG2) << *profile;

    const int length        = profile->num_cols();
    const int num_profiles  = lib_.num_profiles();
    const int center        = lib_.center();
    const int alphabet_size = Alphabet_T::instance().size();

    std::valarray<float> prob(0.0f, num_profiles);  // profile probabilities P(p_k|X_i) at position i
    std::valarray<float> pc(0.0f, alphabet_size);   // pseudocount vector P(a|X_i) at position i
    CountProfile<Alphabet_T>& p = *profile;        // output profile with pseudocounts

    for (int i = 0; i < length; ++i) {
        // calculate profile probabilities P(p_k|X_i)
        for (int k = 0; k < num_profiles; ++k) {
            prob[k] = lib_[k].prior() * pow(2.0, emitter_(lib_[k], p, i));
        }
        prob /= prob.sum();  // normalization

        // calculate pseudocount vector P(a|X_i)
        for(int a = 0; a < alphabet_size; ++a) {
            pc[a] = 0.0f;
            for (int k = 0; k < num_profiles; ++k) {
                pc[a] += prob[k] * pow(2.0, lib_[k][center][a]);
            }
        }
        pc /= pc.sum();  // normalization

        // add pseudocounts to sequence by storing probabilities in output profile
        float tau = pca(p.neff(i));
        for(int a = 0; a < alphabet_size; ++a) {
            p[i][a] = (1.0f - tau) * p[i][a] + tau * pc[a];
        }
    }

    normalize(profile);
    LOG(DEBUG2) << *profile;
}

}  // namespace cs

#endif
