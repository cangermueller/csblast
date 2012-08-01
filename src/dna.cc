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

#include "cs.h"
#include "dna.h"

namespace cs {

// Size of alphabet excluding wildcard character ANY
const size_t Dna::kSize = 4;

// Size of alphabet includding wildcard character ANY
const size_t Dna::kSizeAny = 5;

// Integer code of ANY character
const uint8_t Dna::kAny = 4;

// Integer code of GAP
const uint8_t Dna::kGap = 5;

// Integer code of ENDGAP
const uint8_t Dna::kEndGap = 6;

// For converting from ASCII to DNA integer code
const uint8_t Dna::kCharToInt[] = {
  /*   0 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  16 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  32 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  5,  0,
  /*                                                             -   . */
  /*  48 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  64 */  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
  /*             A   B   C   D   E   F   G   H   I   J   K   L   M   N   O */
  /*  80 */  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*         P   Q   R   S   T   U   V   W   X   Y   Z */
  /*  96 */  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
  /*             a   b   c   d   e   f   g   h   i   j   k   l   m   n   o */
  /* 112 */  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,
  /*         p   q   r   s   t   u   v   w   x   y   z */
  /* 128 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 144 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 160 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 176 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 192 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 208 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 224 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 240 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

// For converting from integer code back to ASCII character
const char Dna::kIntToChar[] = {
  'A', 'C', 'G', 'T', 'N', '-', '-'
};

// For testing if ASCII character is from DNA code
const bool Dna::kValidChar[] = {
  /*   0 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  16 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  32 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,   true,   true,  false,
  /*  48 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  64 */  false,   true,  false,   true,  false,  false,  false,   true,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  80 */  false,  false,  false,  false,   true,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  96 */  false,   true,  false,   true,  false,  false,  false,   true,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 112 */  false,  false,  false,  false,   true,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 128 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 144 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 160 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 176 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 192 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 208 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 224 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 240 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false
};

// Functional groups of DNA alphabet needed for coloring of profile logos
const int Dna::kFuncGroup[] = { 1, 2, 3, 4, 0 };

// Shorthand name for this DNA alphabet
const char Dna::kName[] = "dna";

}  // namespace cs
