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

#ifndef CS_GONNET_MATRIX_H_
#define CS_GONNET_MATRIX_H_

#include "substitution_matrix-inl.h"

namespace cs {

// Gonnet substitution matrix
class GonnetMatrix : public SubstitutionMatrix<AA> {
 public:
  GonnetMatrix();
  virtual ~GonnetMatrix() {}
};

}  // namespace cs

#endif  // CS_GONNET_MATRIX_H_
