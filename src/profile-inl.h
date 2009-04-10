// Copyright 2009, Andreas Biegert

#ifndef SRC_PROFILE_INL_H_
#define SRC_PROFILE_INL_H_

#include "profile.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <limits>
#include <vector>

#include "exception.h"
#include "log.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* Profile<Alphabet>::kClassID = "Profile";

template<class Alphabet>
inline Profile<Alphabet>::Profile()
    : data_(),
      logspace_(false) { }

template<class Alphabet>
inline Profile<Alphabet>::Profile(int num_cols)
    : data_(num_cols, Alphabet::instance().size(), 0.0f),
      logspace_(false) { }

template<class Alphabet>
inline Profile<Alphabet>::Profile(FILE* fin)
    : data_(),
      logspace_(false) {
  read(fin);
}

template<class Alphabet>
inline Profile<Alphabet>::Profile(const Profile& other,
                                  int index,
                                  int length)
    : data_(length, other.alphabet_size(), 0.0f),
      logspace_(other.logspace_) {
  if (index + length > other.num_cols())
    throw Exception("Index=%i and length=%i of sub-profile are out of bounds!",
                    index, length);
  for (int i = 0; i < num_cols(); ++i)
    for (int a = 0; a < alphabet_size(); ++a)
      data_[i][a] = other[i+index][a];
}

template<class Alphabet>
inline void Profile<Alphabet>::readall(FILE* fin,
                                       std::vector< shared_ptr<Profile> >* v) {
  while (!feof(fin)) {
    shared_ptr<Profile> p(new Profile(fin));
    v->push_back(p);

    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
}

template<class Alphabet>
void Profile<Alphabet>::transform_to_logspace() {
  if (!logspace_) {
    for (int i = 0; i < num_cols(); ++i)
      for (int a = 0; a < alphabet_size(); ++a)
        data_[i][a] = log2(data_[i][a]);
    logspace_ = true;
  }
}

template<class Alphabet>
void Profile<Alphabet>::transform_to_linspace() {
  if (logspace_) {
    for (int i = 0; i < num_cols(); ++i)
      for (int a = 0; a < alphabet_size(); ++a)
        data_[i][a] = pow(2.0, data_[i][a]);
    logspace_ = false;
  }
}

template<class Alphabet>
void Profile<Alphabet>::read(FILE* fin) {
  LOG(DEBUG1) << "Reading profile from stream ...";

  // Check if stream actually contains a serialized profile
  char buffer[kBufferSize];
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, class_id()))
    throw Exception("Bad format: profile does not start with '%s'!",
                    class_id());

  read_header(fin);
  read_body(fin);

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void Profile<Alphabet>::read_header(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read num_cols
  int num_cols = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "num_cols")) {
    ptr = buffer;
    num_cols = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'num_cols' record!");
  }
  // Read alphabet_size
  int alphabet_size = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "alphabet_size")) {
    ptr = buffer;
    alphabet_size = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'alphabet_size' record!");
  }
  if (alphabet_size != Alphabet::instance().size())
    throw Exception("Bad format: profile alphabet_size should be %i but is %i!",
                    Alphabet::instance().size(), alphabet_size);
  // Read logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "logspace")) {
    ptr = buffer;
    logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: profile does not contain 'logspace' record!");
  }

  resize(num_cols, alphabet_size);
}

template<class Alphabet>
void Profile<Alphabet>::read_body(FILE* fin) {
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
        data_[i][a] = static_cast<float>(-strtoi_ast(ptr)) / kLogScale;
      else
        data_[i][a] = pow(2.0, static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
    }
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void Profile<Alphabet>::write(FILE* fout) const {
  fprintf(fout, "%s\n", class_id());
  write_header(fout);
  write_body(fout);
}

template<class Alphabet>
void Profile<Alphabet>::write_header(FILE* fout) const {
  fprintf(fout, "num_cols\t%i\n", num_cols());
  fprintf(fout, "alphabet_size\t%i\n", alphabet_size());
  fprintf(fout, "logspace\t%i\n", logspace() ? 1 : 0);
}

template<class Alphabet>
void Profile<Alphabet>::write_body(FILE* fout) const {
  fputc('\t', fout);
  Alphabet::instance().write(fout);

  for (int i = 0; i < num_cols(); ++i) {
    fprintf(fout, "%i", i+1);
    for (int a = 0; a < alphabet_size(); ++a) {
      float log_p = logspace_ ? data_[i][a] : log2(data_[i][a]);
      if (log_p == -std::numeric_limits<float>::infinity())
        fputs("\t*", fout);
      else
        fprintf(fout, "\t%i", -iround(log_p * kLogScale));
    }
    fputc('\n', fout);
  }
  fputs("//\n", fout);
}

template<class Alphabet>
void Profile<Alphabet>::print(std::ostream& out) const {
  out << "\t" << Alphabet::instance() << std::endl;

  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a)
      out << strprintf("\t%6.4f",
                       logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
    out << std::endl;
  }
}

template<class Alphabet>
void Profile<Alphabet>::resize(int num_cols, int alphabet_size) {
  if (num_cols == 0 || alphabet_size == 0)
    throw Exception("Bad profile dimensions: num_cols=%i alphabet_size=%i",
                    num_cols, alphabet_size);
  data_.resize(num_cols, alphabet_size);
}

template<class Alphabet>
inline void reset(Profile<Alphabet>* p) {
  Profile<Alphabet>& profile = *p;
  const int num_cols = profile.num_cols();
  const int alphabet_size = profile.alphabet_size();
  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] = 0.0f;
}

template<class Alphabet>
bool normalize(Profile<Alphabet>* p, float value) {
  Profile<Alphabet>& profile = *p;
  const bool logspace = profile.logspace();
  if (logspace) profile.transform_to_linspace();

  const int num_cols       = profile.num_cols();
  const int alphabet_size  = profile.alphabet_size();
  bool rv = true;

  for (int i = 0; i < num_cols; ++i) {
    float sum = 0.0f;
    for (int a = 0; a < alphabet_size; ++a) sum += profile[i][a];
    if (sum != 0.0f) {
      float fac = value / sum;
      for (int a = 0; a < alphabet_size; ++a) profile[i][a] *= fac;
    } else {
      rv = false;  // couldn't normalize at least one column
    }
  }

  if (logspace) profile.transform_to_logspace();
  return rv;
}

}  // namespace cs

#endif  // SRC_PROFILE_INL_H_
