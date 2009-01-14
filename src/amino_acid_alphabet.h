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
    AminoAcidAlphabet();
    ~AminoAcidAlphabet();

    virtual void unserialize_itoc();

private:
    // Disallow copy and assign
    AminoAcidAlphabet(const AminoAcidAlphabet& );
    AminoAcidAlphabet& operator =(const AminoAcidAlphabet& other);
};

}//cs

#endif
