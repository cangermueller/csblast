/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_TAMURA_NEI_MATRIX_H_
#define CS_TAMURA_NEI_MATRIX_H_

#include "substitution_matrix-inl.h"

namespace cs {

// Substitution matrix according to Tamura & Nei evolutionary model
class TamuraNeiMatrix : public SubstitutionMatrix<Dna> {
 public:
  TamuraNeiMatrix(double time = 0.4,         // evolutionary time
                  const double* f = NULL,    // background freqs (def: 60% AT 40% GC)
                  double alpha_pur = 1.3,    // rate transition purine
                  double alpha_pyr = 1.3,    // rate transition pyrimidine
                  double beta = 1.0) {       // transversion rate

    Vector<double> default_bg_freqs(Dna::kSize, 0.0);
    default_bg_freqs[0] = 0.3;
    default_bg_freqs[1] = 0.2;
    default_bg_freqs[2] = 0.2;
    default_bg_freqs[3] = 0.3;
    if (f == NULL) f = &default_bg_freqs[0];

    for (size_t i = 0; i < Dna::kSize; ++i) {
      for (size_t j = 0; j < Dna::kSize; ++j) {
        const double alphai = (i == 0 || i == 2) ? alpha_pur : alpha_pyr;
        const double pj = (j == 0 || j == 2) ? f[0] + f[2] : f[1] + f[3];
        assert(alphai > 0);
        assert(pj > 0);

        if (i == j) {  // no change
          q_[i][j] = exp(-(alphai + beta) * time)
            + exp(-beta * time) * (1 - exp(-alphai * time)) * (f[j] / pj)
            + (1 - exp(-beta * time)) * f[j];
          q_[i][j] *= f[i]; // convert from subst. prob. to joint prob.

        } else if (abs(static_cast<int>(i) - j) == 2) {  // transition
          q_[i][j] = 0.
            + exp(-beta * time) * (1 - exp(-alphai * time)) * (f[j] / pj)
            + (1 - exp(-beta * time)) * f[j];
            q_[i][j] *= f[i];

        } else { // transversion
          q_[i][j] = 0.
            + 0.
            + (1 - std::exp (-beta * time)) * f[j];
          q_[i][j] *= f[i];
        }
      }
    }

    // Let base class init method do the rest.
    InitFromTargetFreqs();
  }

  virtual ~TamuraNeiMatrix() {}
};

}  // namespace cs

#endif  // CS_BLOSUM_MATRIX_H_
