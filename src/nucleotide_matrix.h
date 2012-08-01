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

#ifndef CS_NUCLEOTIDE_MATRIX_H_
#define CS_NUCLEOTIDE_MATRIX_H_

#include "nucleotide.h"
#include "substitution_matrix-inl.h"

namespace cs {

// BLOSUM family of substitution matrices for  class for substitution matrix
// classes.
class NucleotideMatrix : public SubstitutionMatrix<Nucleotide> {
 public:
  NucleotideMatrix(float match, float mismatch) {
    Init(match, mismatch);
  }
  virtual ~NucleotideMatrix() {}

 private:
  // Initializes all matrix data members.
  virtual void Init(float match, float mismatch);
};

}  // namespace cs

#endif  // CS_NUCLEOTIDE_MATRIX_H_
