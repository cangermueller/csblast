// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_PSEUDOCOUNTS_H_
#define SRC_HMM_PSEUDOCOUNTS_H_

#include <cassert>
#include <cmath>

#include <valarray>

#include "count_profile-inl.h"
#include "mult_emission.h"
#include "hmm-inl.h"
#include "matrix.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from context HMM.
template<class Alphabet>
class HMMPseudocounts : public Pseudocounts<Alphabet> {
 public:
  HMMPseudocounts(const HMM<Alphabet>* hmm,
                  float weight_center,
                  float weight_decay);
  ~HMMPseudocounts() {}

  // Adds context-specific pseudocounts to sequence and stores resulting
  // frequencies in given profile.
  virtual void add_to_sequence(const Sequence<Alphabet>& seq,
                               const Admixture& pca,
                               Profile<Alphabet>* profile) const;
  // Adds context-specific pseudocounts to alignment derived profile.
  virtual void add_to_profile(const Admixture& pca,
                              CountProfile<Alphabet>* p) const;

 private:
  // Profile library with context profiles.
  const HMM<Alphabet>& hmm_;
  // Needed to compute emission probabi
  const MultEmission<Alphabet> emission_;

  DISALLOW_COPY_AND_ASSIGN(HMMPseudocounts);
};  // HMMPseudocounts

}  // namespace cs

#endif  // SRC_HMM_PSEUDOCOUNTS_H_
