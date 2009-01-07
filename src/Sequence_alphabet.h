#ifndef CS_SEQUENCE_ALPHABET_H
#define CS_SEQUENCE_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract interface for alphabet classes whose elements can be represented
// by a sequence of characters (e.g. amino acids or nucleic acids).

#include <cstddef>
#include <vector>
#include <cmath>

class Sequence_alphabet
{
public:
    typedef std::vector<char>::const_iterator const_iterator;

    bool valid(char letter) const;
    size_t size() const;
    int ctoi(char letter) const;
    char itoc(int letter) const;
    int any() const;
    const_iterator begin() const;
    const_iterator end() const;

protected:
    Sequence_alphabet() { }
    ~Sequence_alphabet() { }

    void init();
    virtual std::vector<char> itoc() const = 0;

    std::vector<int> ctoi_;
    std::vector<char> itoc_;

private:
    // Not defined, to prevent copying
    Sequence_alphabet(const Sequence_alphabet& );
    Sequence_alphabet& operator =(const Sequence_alphabet& other);
};


inline bool Sequence_alphabet::valid(char letter) const
{ return ctoi_[letter] != -1; }


inline size_t Sequence_alphabet::size() const
{ return itoc_.size(); }


inline int Sequence_alphabet::ctoi(char letter) const
{ return ctoi_[letter]; }


inline char Sequence_alphabet::itoc(int letter) const
{ return itoc_[letter]; }


inline int Sequence_alphabet::any() const
{ return ctoi_[itoc_[itoc_.size()-1]]; }


inline Sequence_alphabet::const_iterator Sequence_alphabet::begin() const
{ return itoc_.begin(); }


inline Sequence_alphabet::const_iterator Sequence_alphabet::end() const
{ return itoc_.end(); }

#endif
