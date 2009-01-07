#ifndef CS_AMINO_ACID_ALPHABET_H
#define Cs_AMINO_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about an alphabet
// of elements of type T (e.g. amino acids or nucleic acids).

#include "Sequence_alphabet.h"

class Amino_acid_alphabet : public Sequence_alphabet
{
private:
    Amino_acid_alphabet();
    ~Amino_acid_alphabet();

    // Not defined, to prevent copying
    Amino_acid_alphabet(const Amino_acid_alphabet& );
    Amino_acid_alphabet& operator =(const Amino_acid_alphabet& other);

protected:
    virtual std::vector<char> itoc() const;
};

#endif
