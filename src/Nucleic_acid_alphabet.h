#ifndef CS_NUCLEIC_ACID_ALPHABET_H
#define CS_NUCLEIC_ACID_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about nucleic acids

#include "Sequence_alphabet.h"

class Nucleic_acid_alphabet : public Sequence_alphabet
{
public:
    static Nucleic_acid_alphabet* instance();

protected:
    Nucleic_acid_alphabet();
    ~Nucleic_acid_alphabet();

    virtual std::vector<char> itoc() const;

private:
    // Not defined, to prevent copying
    Nucleic_acid_alphabet(const Nucleic_acid_alphabet& );
    Nucleic_acid_alphabet& operator =(const Nucleic_acid_alphabet& other);
};

#endif
