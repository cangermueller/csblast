#ifndef CS_HMM_PSEUDOCOUNTS_H
#define CS_HMM_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of context-specific pseudocounts calculated from context HMM.

#include <cassert>
#include <cmath>

#include <valarray>

#include "count_profile-inl.h"
#include "emitter.h"
#include "forward_backward_algorithm.h"
#include "hmm-inl.h"
#include "log.h"
#include "matrix.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class HMMPseudocounts : public Pseudocounts<Alphabet_T>
{
  public:
    HMMPseudocounts(const HMM<Alphabet_T>* hmm, const EmissionParams& params);
    ~HMMPseudocounts() { }

    // Adds context-specific pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 const Admixture& pca,
                                 Profile<Alphabet_T>* profile) const;
    // Adds context-specific pseudocounts to alignment derived profile.
    virtual void add_to_profile(const Admixture& pca, CountProfile<Alphabet_T>* p) const;

  private:
    // Disallow copy and assign
    HMMPseudocounts(const HMMPseudocounts&);
    void operator=(const HMMPseudocounts&);

    // Profile library with context profiles.
    const HMM<Alphabet_T>& hmm_;
    // Needed to compute emission probabi
    const Emitter<Alphabet_T> emitter_;
};  // HMMPseudocounts



template<class Alphabet_T>
inline HMMPseudocounts<Alphabet_T>::HMMPseudocounts(const HMM<Alphabet_T>* hmm,
                                                    const EmissionParams& params)
        : hmm_(*hmm),
          emitter_(hmm->num_cols(), params)
{
    assert(hmm_.states_logspace());
    assert(!hmm_.transitions_logspace());
}

template<class Alphabet_T>
void HMMPseudocounts<Alphabet_T>::add_to_sequence(const Sequence<Alphabet_T>& seq,
                                                  const Admixture& pca,
                                                  Profile<Alphabet_T>* profile) const
{
    LOG(DEBUG2) << "Adding context-specific, HMM-derived pseudocounts to sequence ...";
    LOG(DEBUG2) << *profile;
    if (seq.length() != profile->num_cols())
        throw Exception("Cannot add context-specific pseudocounts: sequence and profile have different length!");

    const int length        = seq.length();
    const int num_states    = hmm_.num_states();
    const int alphabet_size = Alphabet_T::instance().size();
    const int center        = hmm_.center();
    const float tau         = pca(1.0f);  // number of effective seqs is one

    std::valarray<float> pc(0.0f, alphabet_size);   // pseudocount vector P(a|X_i) at position i
    Profile<Alphabet_T>& p = *profile;              // output profile with pseudocounts

    // calculate posterior state probabilities with forward-backward algorithm
    ForwardBackwardMatrices fbm(length, hmm_.num_states());
    forward_backward_algorithm(hmm_, seq, emitter_, &fbm);

    for (int i = 0; i < length; ++i) {
        // calculate pseudocount vector P(a|X_i)
        for(int a = 0; a < alphabet_size; ++a) {
            pc[a] = 0.0f;
            for (int k = 0; k < num_states; ++k) {
                pc[a] += fbm.f[i][k] * fbm.b[i][k] * pow(2.0, hmm_[k][center][a]);
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
void HMMPseudocounts<Alphabet_T>::add_to_profile(const Admixture& pca, CountProfile<Alphabet_T>* profile) const
{
    assert(!profile->has_counts());
    assert(!profile->logspace());

    LOG(DEBUG2) << "Adding context-specific, HMM-derived pseudocounts to profile ...";
    LOG(DEBUG2) << *profile;

    const int length        = profile->num_cols();
    const int num_states    = hmm_.num_states();
    const int center        = hmm_.center();
    const int alphabet_size = Alphabet_T::instance().size();

    std::valarray<float> pc(0.0f, alphabet_size);   // pseudocount vector P(a|X_i) at position i
    CountProfile<Alphabet_T>& p = *profile;        // output profile with pseudocounts

    // calculate posterior state probabilities with forward-backward algorithm
    ForwardBackwardMatrices fbm(length, hmm_.num_states());
    forward_backward_algorithm(hmm_, p, emitter_, &fbm);

    for (int i = 0; i < length; ++i) {
        // calculate pseudocount vector P(a|X_i)
        for(int a = 0; a < alphabet_size; ++a) {
            pc[a] = 0.0f;
            for (int k = 0; k < num_states; ++k) {
                pc[a] += fbm.f[i][k] * fbm.b[i][k] * pow(2.0, hmm_[k][center][a]);
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

}  // cs

#endif
