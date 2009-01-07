#ifndef CS_AMINO_ACID_ALPHABET_H
#define CS_AMINO_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about amino acids

#include "Sequence_alphabet.h"

class Amino_acid_alphabet : public Sequence_alphabet
{
public:
    static Amino_acid_alphabet* instance();

protected:
    Amino_acid_alphabet();
    ~Amino_acid_alphabet();

    virtual std::vector<char> itoc() const;

private:
    // Not defined, to prevent copying
    Amino_acid_alphabet(const Amino_acid_alphabet& );
    Amino_acid_alphabet& operator =(const Amino_acid_alphabet& other);
};

#endif
