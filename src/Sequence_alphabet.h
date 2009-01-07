#ifndef CS_SEQUENCE_ALPHABET_H
#define Cs_SEQUENCE_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about an alphabet
// of elements of type T (e.g. amino acids or nucleic acids).

#include <cstddef>
#include <vector>
#include <algorithm>
#include <cmath>

class Sequence_alphabet
{
public:
    typedef std::vector<char>::const_iterator const_iterator;

    static Sequence_alphabet* instance()
    {
        // Initialized during first access
        static Sequence_alphabet inst;
        return &inst;
    }

    bool valid(char letter) const;
    size_t size() const;
    int ctoi(char letter) const;
    char itoc(int letter) const;
    int any() const;
    const_iterator begin() const;
    const_iterator end() const;

protected:
    Sequence_alphabet();
    virtual ~Sequence_alphabet();

    // Not defined, to prevent copying
    Sequence_alphabet(const Sequence_alphabet& );
    Sequence_alphabet& operator =(const Sequence_alphabet& other);

    void init();
    virtual std::vector<char> itoc() const;

    std::vector<int> ctoi_;
    std::vector<char> itoc_;
};

#endif
