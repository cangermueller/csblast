// Copyright 2009, Andreas Biegert

#ifndef SRC_CONTEXT_PROFILE_INL_H_
#define SRC_CONTEXT_PROFILE_INL_H_

#include "context_profile.h"

#include <iostream>
#include <vector>

#include "alignment-inl.h"
#include "exception.h"
#include "profile-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* ContextProfile<Alphabet>::kClassID = "ContextProfile";

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile()
    : Profile<Alphabet>(),
      index_(0),
      prior_(0.0f) { }

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(int index, int num_cols)
    : Profile<Alphabet>(num_cols),
      index_(index),
      prior_(0.0f) {
  check();
}

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(int index,
                                                const Profile<Alphabet>& profile)
    : Profile<Alphabet>(profile),
      index_(index),
      prior_(0.0f) {
  check();
}

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(FILE* fin)
    : Profile<Alphabet>(),
      index_(0),
      prior_(0.0f) {
  read(fin);
  check();
}

template<class Alphabet>
void ContextProfile<Alphabet>::check() {
  if (num_cols() % 2 != 1)
    throw Exception("Column number in context profiles must be odd but is %i!",
                    num_cols());
}

template<class Alphabet>
void ContextProfile<Alphabet>::read_header(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read index
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "index")) {
    ptr = buffer;
    index_ = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'index' record!");
  }
  // Read prior
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "prior")) {
    prior_ = atof(buffer + 5);
  } else {
    throw Exception("Bad format: profile does not contain 'prior' record!");
  }

  Profile<Alphabet>::read_header(fin);
}

template<class Alphabet>
void ContextProfile<Alphabet>::write_header(FILE* fout) const {
  fprintf(fout, "index\t\t%i\n", index());
  fprintf(fout, "prior\t\t%-10.8g\n", prior());

  Profile<Alphabet>::write_header(fout);
}

template<class Alphabet>
void ContextProfile<Alphabet>::print(std::ostream& out) const {
  out << "index: " << index() << std::endl;
  out << "prior probability: " << strprintf("%-10.8g", prior()) << std::endl;

  Profile<Alphabet>::print(out);
}

template<class Alphabet>
inline void reset(ContextProfile<Alphabet>* p) {
  ContextProfile<Alphabet>& profile = *p;
  const int num_cols = profile.num_cols();
  const int alphabet_size = profile.alphabet_size();

  profile.set_prior(0.0f);
  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] = 0.0f;
}

}  // namespace cs

#endif  // SRC_CONTEXT_PROFILE_INL_H_
