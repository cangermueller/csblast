// Copyright 2009, Andreas Biegert

#ifndef SRC_CRF_INL_H_
#define SRC_CRF_INL_H_

#include "crf.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "context_weight_state-inl.h"
#include "chain_graph-inl.h"
#include "log.h"
#include "profile-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* CRF<Alphabet>::kClassID = "CRF";

template<class Alphabet>
CRF<Alphabet>::CRF(int num_states, int num_cols)
    : ChainGraph<Alphabet, ContextWeightState>(num_states, num_cols) {}

template<class Alphabet>
CRF<Alphabet>::CRF(FILE* fin)
    : ChainGraph<Alphabet, ContextWeightState>() {
  Read(fin);
}

template<class Alphabet>
CRF<Alphabet>::CRF(
    int num_states,
    int num_cols,
    const StateInitializer<Alphabet, ContextWeightState>& st_init,
    const TransitionInitializer<Alphabet, ContextWeightState>& tr_init)
    : ChainGraph<Alphabet, ContextWeightState>(num_states, num_cols) {
  st_init.Init(*this);
  tr_init.Init(*this);
}

template<class Alphabet>
int CRF<Alphabet>::AddState(const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("CRF contains already %i states!", num_states());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  StatePtr sp(new State(states_.size(), num_states(), profile));

  states_.push_back(sp);
  return states_.size() - 1;
}


template<class Alphabet>
SamplingStateInitializerCRF<Alphabet>::SamplingStateInitializerCRF(
    ProfileVec profiles,
    float sample_rate,
    const Pseudocounts<Alphabet>* pc,
    float pc_admixture)
    : SamplingStateInitializer<Alphabet, ContextWeightState>(profiles,
                                                             sample_rate,
                                                             pc,
                                                             pc_admixture) {}

template<class Alphabet>
LibraryBasedStateInitializerCRF<Alphabet>::LibraryBasedStateInitializerCRF(
    const ProfileLibrary<Alphabet>* lib)
    : LibraryBasedStateInitializer<Alphabet, ContextWeightState>(lib) {}


template<class Alphabet>
void Sparsify(CRF<Alphabet>& crf, float threshold) {
  const bool linspace = !crf.transitions_logspace();
  if (linspace) crf.TransformTransitionsToLogSpace();

  for (int k = 0; k < crf.num_states(); ++k)
    for (int l = 0; l < crf.num_states(); ++l)
      if (crf.test_transition(k,l) && crf(k,l) <= threshold)
        crf.erase_transition(k,l);

  if (linspace) crf.TransformTransitionsToLinSpace();
}

}  // namespace cs

#endif  // SRC_CRF_INL_H_
