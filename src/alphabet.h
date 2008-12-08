#ifndef CS_ALPHABET_H
#define Cs_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about an alphabet
// characters for sequence comparison (e.g. amino acids or nucleic acids).

class Alphabet
{
public:
    static Alphabet* getInstance()
    {
        // Initialized during first access
        static Alphabet inst;
        return &inst;
    }

private:
    Alphabet() {}
    ~Alphabet() {}

    // Not defined, to prevent copying
    Alphabet(const Alphabet& );
    Alphabet& operator =(const Alphabet& other);
};

#endif
