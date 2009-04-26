// Copyright 2009, Andreas Biegert

#ifndef SRC_STATE_INL_H_
#define SRC_STATE_INL_H_

#include "state.h"

#include <cstdio>
#include <cstring>

#include "context_profile-inl.h"
#include "profile-inl.h"
#include "transition.h"

namespace cs {

template<class Alphabet>
const char* State<Alphabet>::kClassID = "State";

template<class Alphabet>
inline State<Alphabet>::State(FILE* fin)
    : ContextProfile<Alphabet>(),
      num_states_(0),
      in_transitions_(0),
      out_transitions_(0) {
  read(fin);
}

template<class Alphabet>
inline State<Alphabet>::State(int index,
                              const Profile<Alphabet>& profile,
                              int num_states)
    : ContextProfile<Alphabet>(index, profile),
      num_states_(num_states),
      in_transitions_(num_states),
      out_transitions_(num_states) {}

template<class Alphabet>
inline State<Alphabet>::State(int index,
                              const ContextProfile<Alphabet>& profile,
                              int num_states)
    : ContextProfile<Alphabet>(profile),
      num_states_(num_states),
      in_transitions_(num_states),
      out_transitions_(num_states) {
  index_ = index;
}

template<class Alphabet>
inline float State<Alphabet>::to(int k) const {
  return out_transitions_.test(k) ? out_transitions_.get(k).probability : 0.0f;
}

template<class Alphabet>
inline float State<Alphabet>::from(int k) const {
  return in_transitions_.test(k) ? in_transitions_.get(k).probability : 0.0f;
}

template<class Alphabet>
void State<Alphabet>::clear_transitions() {
  in_transitions_.clear();
  out_transitions_.clear();
}

template<class Alphabet>
void State<Alphabet>::resize(int num_states) {
  clear_transitions();
  in_transitions_.resize(num_states);
  out_transitions_.resize(num_states);
}

template<class Alphabet>
void State<Alphabet>::read_header(FILE* fin) {
  ContextProfile<Alphabet>::read_header(fin);

  // Read HMM size
  char buffer[kBufferSize];
  const char* ptr = buffer;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states_ = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'NSTATES' record!");
  }

  in_transitions_.resize(num_states_);
  out_transitions_.resize(num_states_);
}

template<class Alphabet>
void State<Alphabet>::write_header(FILE* fout) const {
  ContextProfile<Alphabet>::write_header(fout);
  fprintf(fout, "NSTATES\t%i\n", num_states_);
}

}  // namespace cs

#endif  // SRC_STATE_INL_H_
