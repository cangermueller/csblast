// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_STATE_INL_H_
#define SRC_HMM_STATE_INL_H_

#include "hmm_state.h"

#include <cstdio>
#include <cstring>

#include "context_profile-inl.h"
#include "profile-inl.h"
#include "transition.h"

namespace cs {

template<class Alphabet>
const char* HMMState<Alphabet>::kClassID = "HMMState";

template<class Alphabet>
HMMState<Alphabet>::HMMState(FILE* fin)
    : ContextProfile<Alphabet>() {
  Read(fin);
}

template<class Alphabet>
HMMState<Alphabet>::HMMState(int index,
                             int num_states,
                             const Profile<Alphabet>& profile)
    : ContextProfile<Alphabet>(index, profile),
      in_transitions_(num_states),
      out_transitions_(num_states) {}

template<class Alphabet>
HMMState<Alphabet>::HMMState(int index,
                             int num_states,
                             const ContextProfile<Alphabet>& profile)
    : ContextProfile<Alphabet>(profile),
      in_transitions_(num_states),
      out_transitions_(num_states) {
  index_ = index;
}

template<class Alphabet>
void HMMState<Alphabet>::ClearTransitions() {
  in_transitions_.clear();
  out_transitions_.clear();
}

template<class Alphabet>
void HMMState<Alphabet>::ReadHeader(FILE* fin) {
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
void HMMState<Alphabet>::WriteHeader(FILE* fout) const {
  ContextProfile<Alphabet>::WriteHeader(fout);
  fprintf(fout, "NSTATES\t%i\n", static_cast<int>(in_transitions_.size()));
}

}  // namespace cs

#endif  // SRC_HMM_STATE_INL_H_
