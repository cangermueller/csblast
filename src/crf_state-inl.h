// Copyright 2009, Andreas Biegert

#ifndef SRC_CRF_STATE_INL_H_
#define SRC_CRF_STATE_INL_H_

#include "crf_state.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "log.h"
#include "profile-inl.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
CRFState<Alphabet>::CRFState()
    : index_(0) {}

template<class Alphabet>
CRFState<Alphabet>::CRFState(int index,
                            int num_states,
                            const Profile<Alphabet>& profile)
    : index_(index),
      in_transitions_(num_states),
      out_transitions_(num_states) {
  Init(profile);
}

template<class Alphabet>
void CRFState<Alphabet>::Read(FILE* fin) {
  // Check if stream actually contains a serialized profile
  char buffer[kBufferSize];
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, "CRFState"))
    throw Exception("Bad format: profile does not start with 'CRFState'!");

  ReadHeader(fin);
  ReadBody(fin);

  //  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void CRFState<Alphabet>::ReadHeader(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read index
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "INDEX")) {
    ptr = buffer;
    index_ = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'INDEX' record!");
  }
  // Read num_states
  int num_states = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'NSTATES' record!");
  }
  // Read num_cols
  int num_cols = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NCOLS")) {
    ptr = buffer;
    num_cols = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'NCOLS' record!");
  }
  // Read alphabet_size
  int alphabet_size = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ALPH")) {
    ptr = buffer;
    alphabet_size = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'ALPH' record!");
  }
  if (alphabet_size != Alphabet::instance().size())
    throw Exception("Bad format: profile alphabet_size should be %i but is %i!",
                    Alphabet::instance().size(), alphabet_size);

  Resize(num_cols, alphabet_size);
  in_transitions_.resize(num_states);
  out_transitions_.resize(num_states);
}

template<class Alphabet>
void CRFState<Alphabet>::ReadBody(FILE* fin) {
  const int alph_size = alphabet_size();
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int i = 0;

  fgetline(buffer, kBufferSize, fin);  // skip alphabet description line
  while (fgetline(buffer, kBufferSize, fin)
         && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    i = strtoi(ptr) - 1;
    for (int a = 0; a < alph_size; ++a) {
      if (logspace_)
        weights_[i][a] = static_cast<float>(-strtoi_ast(ptr)) / kLogScale;
      else
        weights_[i][a] = fast_pow2(static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
    }
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void CRFState<Alphabet>::Write(FILE* fout) const {
  fprintf(fout, "%s\n", class_id());
  WriteHeader(fout);
  WriteBody(fout);
}

template<class Alphabet>
void CRFState<Alphabet>::WriteHeader(FILE* fout) const {
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ALPH\t%i\n", alphabet_size());
  fprintf(fout, "LOG\t%i\n", logspace() ? 1 : 0);
}

template<class Alphabet>
void CRFState<Alphabet>::WriteBody(FILE* fout) const {
  fputs("WGHT\t", fout);
  Alphabet::instance().Write(fout);
  fputc('\n', fout);

  for (int i = 0; i < num_cols(); ++i) {
    fprintf(fout, "%i", i+1);
    for (int a = 0; a < alphabet_size(); ++a) {
      fprintf(fout, "\t%+i", -iround(weights_[i][a] * kLogScale));
    }
    fputc('\n', fout);
  }

  fputs("WGHT", fout);
  for (int a = 0; a < alphabet_size(); ++a)
    fprintf(fout, "\t%+i", -iround(pc_[a] * kLogScale));
  fputs("//\n", fout);
}

template<class Alphabet>
void CRFState<Alphabet>::Print(std::ostream& out) const {
  out << "\t" << Alphabet::instance() << std::endl;

  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a)
      out << strprintf("\t%6.2f", weights_[i][a]);
    out << std::endl;
  }
}

template<class Alphabet>
void CRFState<Alphabet>::Resize(int num_cols, int alphabet_size) {
  if (num_cols == 0 || alphabet_size == 0)
    throw Exception("Bad profile dimensions: num_cols=%i alphabet_size=%i",
                    num_cols, alphabet_size);
  weights_.resize(num_cols, alphabet_size);
}

template<class Alphabet>
inline void Reset(CRFState<Alphabet>* p) {
  CRFState<Alphabet>& profile = *p;
  const int num_cols = profile.num_cols();
  const int alphabet_size = profile.alphabet_size();
  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] = 0.0f;
}

}  // namespace cs

#endif  // SRC_PROFILE_INL_H_
