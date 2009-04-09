// Copyright 2009, Andreas Biegert

#ifndef SRC_STATE_INL_H_
#define SRC_STATE_INL_H_

#include "state.h"

#include <string>

namespace cs {

template<class Alphabet>
const char* State<Alphabet>::kClassID = "State";

template<class Alphabet>
inline State<Alphabet>::State(std::istream& in)
    : ContextProfile<Alphabet>(),
      num_states_(0),
      in_transitions_(0),
      out_transitions_(0) {
  read(in);
}

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
      out_transitions_(num_states) { }

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
void State<Alphabet>::read_header(std::istream& in) {
  ContextProfile<Alphabet>::read_header(in);

  // Read HMM size
  std::string tmp;
  if (getline(in, tmp) && tmp.find("num_states") != std::string::npos)
    num_states_ = atoi(tmp.c_str() + 10);

  in_transitions_.resize(num_states_);
  out_transitions_.resize(num_states_);
}

template<class Alphabet>
void State<Alphabet>::read_header(FILE* fin) {
  ContextProfile<Alphabet>::read_header(fin);

  // Read HMM size
  char buffer[kBufferSize];
  const char* ptr = buffer;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "num_states")) {
    ptr = buffer;
    num_states_ = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'num_states' record!");
  }

  in_transitions_.resize(num_states_);
  out_transitions_.resize(num_states_);
}

template<class Alphabet>
void State<Alphabet>::write_header(std::ostream& out) const {
  ContextProfile<Alphabet>::write_header(out);
  out << "num_states\t" << num_states_ << std::endl;
}

template<class Alphabet>
void State<Alphabet>::write_header(FILE* fout) const {
  ContextProfile<Alphabet>::write_header(fout);

  fprintf(fout, "num_states\t%i\n", num_states_);
}

}  // namespace cs

#endif  // SRC_STATE_INL_H_
