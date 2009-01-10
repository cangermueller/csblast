#ifndef CS_NUCLEIC_ACID_ALPHABET_H
#define CS_NUCLEIC_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about nucleic acids

#include "sequence_alphabet.h"
#include <iostream>

namespace cs
{

class NucleicAcidAlphabet : public SequenceAlphabet
{
public:
    static NucleicAcidAlphabet* instance();

protected:
    NucleicAcidAlphabet();
    ~NucleicAcidAlphabet();

    virtual void init_itoc();

private:
    // Not defined, to prevent copying
    NucleicAcidAlphabet(const NucleicAcidAlphabet& );
    NucleicAcidAlphabet& operator =(const NucleicAcidAlphabet& other);
};

}//cs

#endif
