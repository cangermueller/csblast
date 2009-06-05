// Copyright 2009, Andreas Biegert

#ifndef SRC_CONTEXT_PROFILE_STATE_INL_H_
#define SRC_CONTEXT_PROFILE_STATE_INL_H_

#include "context_profile_state.h"

#include <cstdio>
#include <cstring>

#include "context_profile-inl.h"
#include "profile-inl.h"
#include "transition.h"

namespace cs {

template<class Alphabet>
const char* ContextProfileState<Alphabet>::kClassID = "ContextProfileState";

template<class Alphabet>
ContextProfileState<Alphabet>::ContextProfileState(FILE* fin)
    : ContextProfile<Alphabet>() {
  Read(fin);
}

template<class Alphabet>
ContextProfileState<Alphabet>::ContextProfileState(int index,
                             int num_states,
                             const Profile<Alphabet>& profile)
    : ContextProfile<Alphabet>(index, profile),
      in_transitions_(num_states),
      out_transitions_(num_states) {}

template<class Alphabet>
void ContextProfileState<Alphabet>::ClearTransitions() {
  in_transitions_.clear();
  out_transitions_.clear();
}

template<class Alphabet>
void ContextProfileState<Alphabet>::ReadHeader(FILE* fin) {
  ContextProfile<Alphabet>::ReadHeader(fin);

  // Read HMM size
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int num_states = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'NSTATES' record!");
  }

  in_transitions_.resize(num_states);
  out_transitions_.resize(num_states);
}

template<class Alphabet>
void ContextProfileState<Alphabet>::WriteHeader(FILE* fout) const {
  ContextProfile<Alphabet>::WriteHeader(fout);
  fprintf(fout, "NSTATES\t%i\n", static_cast<int>(in_transitions_.size()));
}

}  // namespace cs

#endif  // SRC_CONTEXT_PROFILE_STATE_INL_H_
