// Copyright 2009, Andreas Biegert

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
