// Copyright 2009, Andreas Biegert

#ifndef SRC_AMINO_ACID_H_
#define SRC_AMINO_ACID_H_

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

#endif  // SRC_AMINO_ACID_H_
