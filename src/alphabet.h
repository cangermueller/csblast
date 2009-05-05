// Copyright 2009, Andreas Biegert

#ifndef SRC_ALPHABET_H_
#define SRC_ALPHABET_H_

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <vector>

#include "globals.h"

namespace cs {

// Abstract base class for alphabet classes whose elements can be represented
// by a sequence of characters (e.g. amino acids or nucleic acids).
class Alphabet {
 public:
  typedef std::vector<char>::const_iterator const_iterator;

  // Returns the number of letters in the alphabet (incl. ANY)
  int size() const { return size_; }
  // Returns the integer representation of the given character.
  int ctoi(char letter) const { return ctoi_[toupper(letter)]; }
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
  // Returns true if the character belongs to the alphabet.
  bool valid(char letter, bool allow_gap = false) const {
    return ctoi_[toupper(letter)] != kInvalidChar
      && (allow_gap || (letter != '-' && letter != '.'));
  }
  // Returns a const iterator to the first character in the alphabet.
  const_iterator begin() const { return itoc_.begin(); }
  // Returns a const iterator just past the end of last distinct character in
  // the alphabet.
  const_iterator end() const { return itoc_.begin() + size(); }

  // Print alphabet characters delimited by tabs.
  friend std::ostream& operator<< (std::ostream& out, const Alphabet& alph) {
    out << *alph.begin();
    for (Alphabet::const_iterator a = alph.begin() + 1; a != alph.end(); ++a)
      out << '\t' << *a;
    return out;
  }
  // Write alphabet characters delimited by tabs.
  void write(FILE* fout) const;

 protected:
  // Denotes invalid characters in ctoi array
  static const int kInvalidChar = -1;

  // Constructor to be used by derived classes to setup alphabet.
  Alphabet(int size, char any);
  virtual ~Alphabet() {}

  // Initializes ctoi and itoc conversion arrays.
  void Init();
  // Allows derived classes to set additional ctoi conversions.
  void set_ctoi(char letter, int i) { ctoi_[toupper(letter)] = i; }

 private:
  // Gets ctoi conversion array from derived class.
  virtual const char* get_itoc() const = 0;

  // Size of alphabet (incl. ANY character).
  const int size_;
  // Any character of alphabet in character representation.
  const char any_;
  // Conversion array from character to iteger representation.
  std::vector<int> ctoi_;
  // Conversion array from integer to character representation
  // (incl. ANY, GAP, and ENDGAP).
  std::vector<char> itoc_;

  DISALLOW_COPY_AND_ASSIGN(Alphabet);
};

}  // namespace cs

#endif  // SRC_ALPHABET_H_
