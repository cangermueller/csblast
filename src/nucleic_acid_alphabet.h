#ifndef CS_NUCLEIC_ACID_ALPHABET_H
#define CS_NUCLEIC_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about nucleic acids

#include "sequence_alphabet.h"

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
    // Disallow copy and assign
    NucleicAcidAlphabet(const NucleicAcidAlphabet& );
    NucleicAcidAlphabet& operator =(const NucleicAcidAlphabet& other);
};

}//cs

#endif
