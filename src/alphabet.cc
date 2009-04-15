// Copyright 2009, Andreas Biegert

#include "alphabet.h"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "exception.h"

namespace cs {

Alphabet::Alphabet(int size, char any)
    : size_(size),
      any_(any),
      ctoi_(static_cast<int>(pow(2, 8*sizeof(char))), kInvalidChar),
      itoc_(size + 3, '\0') {}

void Alphabet::init() {
  const char* itoc = get_itoc();
  if (static_cast<int>(strlen(itoc)) != size_)
    throw Exception("Alphabet error: itoc array as length %i but should have %i!",
                    strlen(itoc), size_);

  // Setup itoc vector
  copy(itoc, itoc + strlen(itoc), itoc_.begin());
  itoc_[size_]   = any_;   // ANY
  itoc_[size_ + 1] = '-';  // GAP
  itoc_[size_ + 2] = '-';  // ENDGAP

  // Setup ctoi vector
  for (int i = 0; i < size_; ++i) ctoi_[toupper(itoc_[i])] = i;
  ctoi_[toupper(any_)] = size_;      // ANY
  ctoi_['-']           = size_ + 1;  // MATCH GAP
  ctoi_['.']           = size_ + 1;  // INSERT GAP
}

void Alphabet::write(FILE* fout) const {
  fputc(itoc_[0], fout);
  for (int a = 1; a < size_; ++a) {
    fputc('\t', fout);
    fputc(itoc_[a], fout);
  }
}

}  // namespace cs
