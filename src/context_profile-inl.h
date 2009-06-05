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
      prior_(0.0) {}

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(int index, int num_cols)
    : Profile<Alphabet>(num_cols),
      index_(index),
      prior_(0.0) {
  check();
}

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(int index,
                                                const Profile<Alphabet>& profile)
    : Profile<Alphabet>(profile),
      index_(index),
      prior_(0.0) {
  check();
}

template<class Alphabet>
inline ContextProfile<Alphabet>::ContextProfile(FILE* fin)
    : Profile<Alphabet>(),
      index_(0),
      prior_(0.0) {
  Read(fin);
  check();
}

template<class Alphabet>
void ContextProfile<Alphabet>::check() {
  if (num_cols() % 2 != 1)
    throw Exception("Column number in context profiles must be odd but is %i!",
                    num_cols());
}

template<class Alphabet>
void ContextProfile<Alphabet>::ReadHeader(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read index
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "INDEX")) {
    ptr = buffer;
    index_ = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'INDEX' record!");
  }
  // Read prior
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "PRIOR")) {
    prior_ = atof(buffer + 5);
  } else {
    throw Exception("Bad format: profile does not contain 'PRIOR' record!");
  }

  Profile<Alphabet>::ReadHeader(fin);
}

template<class Alphabet>
void ContextProfile<Alphabet>::WriteHeader(FILE* fout) const {
  fprintf(fout, "INDEX\t%i\n", index());
  fprintf(fout, "PRIOR\t%-10.8g\n", prior());

  Profile<Alphabet>::WriteHeader(fout);
}

template<class Alphabet>
void ContextProfile<Alphabet>::Print(std::ostream& out) const {
  out << "index: " << index() << std::endl;
  out << "prior probability: " << strprintf("%-10.8g", prior()) << std::endl;

  Profile<Alphabet>::Print(out);
}

template<class Alphabet>
inline void Reset(ContextProfile<Alphabet>* p) {
  ContextProfile<Alphabet>& profile = *p;
  const int num_cols = profile.num_cols();
  const int alphabet_size = profile.alphabet_size();

  profile.set_prior(0.0);
  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] = 0.0f;
}

}  // namespace cs

#endif  // SRC_CONTEXT_PROFILE_INL_H_
