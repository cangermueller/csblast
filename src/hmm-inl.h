// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_INL_H_
#define SRC_HMM_INL_H_

#include "hmm.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "co_emission-inl.h"
#include "context_profile_state-inl.h"
#include "chain_graph-inl.h"
#include "count_profile-inl.h"
#include "log.h"
#include "profile-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* HMM<Alphabet>::kClassID = "HMM";

template<class Alphabet>
HMM<Alphabet>::HMM(int num_states, int num_cols)
    : ChainGraph<Alphabet, ContextProfileState>(num_states, num_cols),
      states_logspace_(false)  {}

template<class Alphabet>
HMM<Alphabet>::HMM(FILE* fin)
    : ChainGraph<Alphabet, ContextProfileState>(),
      states_logspace_(false) {
  Read(fin);
}

template<class Alphabet>
HMM<Alphabet>::HMM(
    int num_states,
    int num_cols,
    const StateInitializer<Alphabet, ContextProfileState>& st_init,
    const TransitionInitializer<Alphabet, ContextProfileState>& tr_init)
    : ChainGraph<Alphabet, ContextProfileState>(num_states, num_cols),
      states_logspace_(false) {
  st_init.Init(*this);
  tr_init.Init(*this);
}

template<class Alphabet>
int HMM<Alphabet>::AddState(const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("HMM contains already %i states!", num_states());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  StatePtr sp(new ContextProfileState<Alphabet>(states_.size(),
                                                num_states(),
                                                profile));
  sp->set_prior(1.0 / num_states());
  states_.push_back(sp);

  return states_.size() - 1;
}

template<class Alphabet>
inline void HMM<Alphabet>::TransformStatesToLogSpace() {
  if (!states_logspace()) {
    for (StateIter si = states_begin(); si != states_end(); ++si)
      (*si)->TransformToLogSpace();
    states_logspace_ = true;
  }
}

template<class Alphabet>
inline void HMM<Alphabet>::TransformStatesToLinSpace() {
  if (states_logspace()) {
    for (StateIter si = states_begin(); si != states_end(); ++si)
      (*si)->TransformToLinSpace();
    states_logspace_ = false;
  }
}

template<class Alphabet>
void HMM<Alphabet>::ReadHeader(FILE* fin) {
  ChainGraph<Alphabet, ContextProfileState>::ReadHeader(fin);

  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read states logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "STLOG")) {
    ptr = buffer;
    states_logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: HMM does not contain 'STLOG' record!");
  }
}

template<class Alphabet>
void HMM<Alphabet>::WriteHeader(FILE* fout) const {
  ChainGraph<Alphabet, ContextProfileState>::WriteHeader(fout);
  fprintf(fout, "STLOG\t%i\n", states_logspace() ? 1 : 0);
}


template<class Alphabet>
SamplingStateInitializerHMM<Alphabet>::SamplingStateInitializerHMM(
    ProfileVec profiles,
    float sample_rate,
    const Pseudocounts<Alphabet>* pc,
    float pc_admixture)
    : SamplingStateInitializer<Alphabet, ContextProfileState>(profiles,
                                                              sample_rate,
                                                              pc,
                                                              pc_admixture) {}

template<class Alphabet>
LibraryStateInitializerHMM<Alphabet>::LibraryStateInitializerHMM(
    const ProfileLibrary<Alphabet>* lib)
    : LibraryStateInitializer<Alphabet, ContextProfileState>(lib) {}

template<class Alphabet>
CoEmissionTransitionInitializerHMM<Alphabet>::CoEmissionTransitionInitializerHMM(
    const SubstitutionMatrix<Alphabet>* sm,
    float thresh)
    : CoEmissionTransitionInitializer<Alphabet, ContextProfileState>(sm, thresh) {}

}  // namespace cs

#endif  // SRC_HMM_INL_H_
