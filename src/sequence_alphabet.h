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

namespace cs
{

class SequenceAlphabet
{
public:
    typedef std::vector<char>::const_iterator const_iterator;

    bool valid(char letter) const;
    int size() const;
    int ctoi(char letter) const;
    char itoc(int letter) const;
    int any() const;
    const_iterator begin() const;
    const_iterator end() const;

protected:
    SequenceAlphabet();
    ~SequenceAlphabet();

    // Initializes ctoi and itoc conversion arrays.
    void init();
    // Template method that initializes itoc conversion array.
    // Note: The last element of itoc has to specify the "any" character.
    virtual void init_itoc() = 0;

    static const int kInvalidChar = -1;
    std::vector<int> ctoi_;
    std::vector<char> itoc_;

private:
    // Not defined, to prevent copying
    SequenceAlphabet(const SequenceAlphabet& );
    SequenceAlphabet& operator =(const SequenceAlphabet& other);
};



inline bool SequenceAlphabet::valid(char letter) const
{ return ctoi_[letter] != kInvalidChar; }

inline int SequenceAlphabet::size() const
{ return itoc_.size(); }

inline int SequenceAlphabet::ctoi(char letter) const
{ return ctoi_[static_cast<int>(letter)]; }

inline char SequenceAlphabet::itoc(int letter) const
{ return itoc_[letter]; }

inline int SequenceAlphabet::any() const
{ return ctoi_[itoc_[itoc_.size()-1]]; }

inline SequenceAlphabet::const_iterator SequenceAlphabet::begin() const
{ return itoc_.begin(); }

inline SequenceAlphabet::const_iterator SequenceAlphabet::end() const
{ return itoc_.end(); }

}//cs

#endif
