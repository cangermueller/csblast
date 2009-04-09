// Copyright 2009, Andreas Biegert

#include "nucleotide_matrix.h"

namespace cs {

void NucleotideMatrix::init(float match, float mismatch) {
  // Set background vector
  f_.assign(size_, 1.0f / size_);

  // Fill substitution matrix
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b < size_; ++b)
      s_[a][b] = a==b ? match : mismatch;

  init_from_scores_and_background_frequencies();
}

}  // namespace cs
