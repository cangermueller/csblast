#ifndef CS_SEQUENCE_ALPHABET_H
#define CS_SEQUENCE_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract interface for alphabet classes whose elements can be represented
// by a sequence of characters (e.g. amino acids or nucleic acids).

#include <vector>

namespace cs
{

class SequenceAlphabet
{
public:
    typedef std::vector<char>::const_iterator const_iterator;

    SequenceAlphabet() {}
    virtual ~SequenceAlphabet() {}

    // Returns true if the character belongs to the alphabet.
    bool valid(char letter) const { return ctoi_[letter] != kInvalidChar; }
    // Returns the number of letters in the alphabet (incl. ANY)
    int size() const { return itoc_.size(); }
    // Returns the integer representation of the given character.
    int ctoi(char letter) const { return ctoi_[static_cast<int>(letter)]; }
    // Returns the character representation of the given integer.
    char itoc(int letter) const { return itoc_[letter]; }
    int any() const { return ctoi_[itoc_[itoc_.size()-1]]; }
    const_iterator begin() const { return itoc_.begin(); }
    const_iterator end() const { return itoc_.end(); }

protected:
    // Denotes invalid characters in ctoi array
    static const int kInvalidChar = -1;

    // Initializes ctoi and itoc conversion arrays.
    void init();
    // Initializes itoc conversion array. Should be implemented by derived classes.
    // Note: Last element in itoc_ should be the ANY character.
    virtual void init_itoc() = 0;

    // Conversion array from character to iteger representation.
    std::vector<int> ctoi_;
    // Conversion array from integer to character representation. Last element is
    // ANY character.
    std::vector<char> itoc_;

private:
    // Disallow copy and assign
    SequenceAlphabet(const SequenceAlphabet& other);
    SequenceAlphabet& operator =(const SequenceAlphabet& other);
};

}//cs

#endif
