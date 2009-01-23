#ifndef CS_NUCLEOTIDE_ALPHABET_H
#define CS_NUCLEOTIDE_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about nucleic acids

#include "sequence_alphabet.h"

namespace cs
{

class NucleotideAlphabet : public SequenceAlphabet
{
  public:
    static NucleotideAlphabet* instance();

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return nucleotides_; }

  private:
    // IUPAC nucleotide code
    static const char nucleotides_[];

    NucleotideAlphabet();
    ~NucleotideAlphabet() {}
    // Disallow copy and assign
    NucleotideAlphabet(const NucleotideAlphabet& );
    NucleotideAlphabet& operator =(const NucleotideAlphabet& other);
};

}  // cs

#endif
