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

#ifndef CS_AMINO_ACID_H_
#define CS_AMINO_ACID_H_

#include "alphabet.h"
#include "globals.h"

namespace cs {

// Singleton class that encapsulates meta information about amino acids.
class AminoAcid : public Alphabet {
 public:

  // Returns an instance of the amino acid alphabet singleton.
  static const AminoAcid& instance() {
    static AminoAcid inst;
    return inst;
  }

 protected:
  // Gets itoc conversion array from derived class.
  virtual const char* get_itoc() const { return amino_acids_; }

 private:
  // IUPAC amino acid code
  static const char amino_acids_[];

  AminoAcid() : Alphabet(20, 'X') {
    Init();
    // Conversion of non-standard amino acids to integers
    set_ctoi('O', any());
    set_ctoi('J', any());
    set_ctoi('U', ctoi('C'));
    set_ctoi('B', ctoi('D'));
    set_ctoi('Z', ctoi('E'));
  }
  virtual ~AminoAcid() {}

  DISALLOW_COPY_AND_ASSIGN(AminoAcid);
};

}  // namespace cs

#endif  // CS_AMINO_ACID_H_
