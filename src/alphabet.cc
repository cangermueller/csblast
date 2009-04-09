// Copyright 2009, Andreas Biegert

#include "alphabet.h"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>

#include "exception.h"

namespace cs {

Alphabet::Alphabet(int size, char any)
        : size_(size),
          any_(any),
          ctoi_(static_cast<int>(pow(2, 8*sizeof(char))), INVALID_CHAR),
          itoc_(size + 3, '\0') { }

void Alphabet::init() {
  const char* itoc = get_itoc();
  if (static_cast<int>(strlen(itoc)) != size_)
    throw Exception("Int2Char array has length %i but should have length %i!",
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

}  // namespace cs
