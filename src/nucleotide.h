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

#ifndef CS_NUCLEOTIDE_H_
#define CS_NUCLEOTIDE_H_

#include "alphabet.h"
#include "globals.h"

namespace cs {

// Singleton class that encapsulates meta information about nucleic acids.
class Nucleotide : public Alphabet {
 public:

  // Returns an instance of the nucleotide alphabet singleton.
  static const Nucleotide& instance() {
    static Nucleotide inst;
    return inst;
  }

 protected:
  // Gets ctoi conversion array from derived class.
  virtual const char* get_itoc() const { return nucleotides_; }

 private:
  // IUPAC nucleotide code
  static const char nucleotides_[];

  Nucleotide() : Alphabet(4, 'N') { Init(); }
  virtual ~Nucleotide() {}

  DISALLOW_COPY_AND_ASSIGN(Nucleotide);
};

}  // namespace cs

#endif  // CS_NUCLEOTIDE_H_
