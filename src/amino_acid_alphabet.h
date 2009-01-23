#ifndef CS_AMINO_ACID_ALPHABET_H
#define CS_AMINO_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about amino acids

#include "sequence_alphabet.h"

namespace cs
{

class AminoAcidAlphabet : public SequenceAlphabet
{
  public:
    static AminoAcidAlphabet* instance();

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return amino_acids_; }

  private:
    // IUPAC amino acid code
    static const char amino_acids_[];

    AminoAcidAlphabet();
    ~AminoAcidAlphabet() {}
    // Disallow copy and assign
    AminoAcidAlphabet(const AminoAcidAlphabet& );
    AminoAcidAlphabet& operator =(const AminoAcidAlphabet& other);
};

}  // cs

#endif
