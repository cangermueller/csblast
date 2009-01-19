#ifndef CS_SEQUENCE_ALPHABET_H
#define CS_SEQUENCE_ALPHABET_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for alphabet classes whose elements can be represented
// by a sequence of characters (e.g. amino acids or nucleic acids).

#include <vector>

namespace cs
{

class SequenceAlphabet
{
  public:
    // Returns the number of letters in the alphabet (incl. ANY)
    int size() const { return size_; }
    // Returns the integer representation of the given character.
    int ctoi(char letter) const { return ctoi_[static_cast<int>(letter)]; }
    // Returns the character representation of the given integer.
    char itoc(int letter) const { return itoc_[letter]; }
    // Returns the integer representation of the any character.
    int any() const { return size_; }
    // Returns integer representation of gap.
    int gap() const { return size_ + 1; }
    // Returns integer representation of endgap.
    int endgap() const { return size_ + 2; }
    // Returns true if provided integer represents ANY
    bool any(int letter) const { return letter == size_; }
    // Returns true if provided integer represents a real symbol (letter < ANY)
    bool less_any(int letter) const { return letter < size_; }
    // Returns true if provided integer represents GAP
    bool gap(int letter) const { return letter == size_ + 1; }
    // Returns true if provided integer represents ENDGAP
    bool endgap(int letter) const { return letter == size_ + 2; }
    // Returns the any character.
    char any_chr() const { return any_; }
    // Returns the gap character.
    char gap_chr() const { return gap_; }
    // Returns true if the character belongs to the alphabet.
    bool valid(char letter, bool allow_gap = false) const
    { return ctoi_[letter] != kInvalidChar && (allow_gap || letter != gap_); }

  protected:
    // Denotes invalid characters in ctoi array
    static const int kInvalidChar = -1;

    // Constructor to be used by derived classes to setup alphabet.
    SequenceAlphabet(int size, char any, char gap = '-');
    ~SequenceAlphabet() {}

    // Initializes ctoi and itoc conversion arrays.
    void init();
    // Allows derived classes to set additional ctoi conversions.
    void set_ctoi(char c, int i) { ctoi_[c] = i; }

  private:
    // Disallow copy and assign
    SequenceAlphabet(const SequenceAlphabet& other);
    SequenceAlphabet& operator =(const SequenceAlphabet& other);

    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const = 0;

    // Size of alphabet (incl. ANY character).
    const int size_;
    // Any character of alphabet in character representation.
    const char any_;
    // Gap character of alphabet in character representation.
    const char gap_;
    // Conversion array from character to iteger representation.
    std::vector<int> ctoi_;
    // Conversion array from integer to character representation (incl. ANY, GAP, and ENDGAP).
    std::vector<char> itoc_;
};

}  // cs

#endif
