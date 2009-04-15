// Copyright 2009, Andreas Biegert

#ifndef SRC_SEQUENCE_H_
#define SRC_SEQUENCE_H_

#include <cctype>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "globals.h"
#include "shared_ptr.h"

namespace cs {

// A container class representing a sequence consisting of letters over a
// sequence alphabet.
template<class Alphabet>
class Sequence {
 public:
  typedef char* iterator;
  typedef const char* const_iterator;

  // Constructs sequence with specified length.
  explicit Sequence(int length);
  // Constructs sequence from serialized sequence in FASTA format read from
  // file stream.
  explicit Sequence(FILE* fin);
  // Constructs sequence with given header and sequence string of characters.
  Sequence(const std::string& header, const std::string& sequence);

  ~Sequence() {}

  // Reads all available sequences from the input stream and returns them in a
  // vector.
  static void readall(FILE* fin, std::vector< shared_ptr<Sequence> >* v);

  // Accessors for integer at position i of the sequence.
  char& operator[](int i) { return seq_[i]; }
  const char& operator[](int i) const { return seq_[i]; }
  char& at(int i) { return seq_[i]; }
  // Returns the character at position i of the sequence.
  char chr(int i) const { return Alphabet::instance().itoc(seq_[i]); }
  // Returns the sequence length.
  int length() const { return seq_.size(); }
  // Returns the header string of the sequence.
  std::string header() const { return header_; }
  // Sets the header to given string.
  void set_header(const std::string& header) { header_ = header; }
  // Returns a const iterator to the first integer element of the sequence.
  const_iterator begin() const { return &seq_[0]; }
  // Returns a const iterator just past the end of the sequence.
  const_iterator end() const { return begin() + length(); }
  // Returns an iterator to the first integer element of the sequence.
  iterator begin() { return &seq_[0]; }
  // Returns an iterator just past the end of the sequence.
  iterator end() { return begin() + length(); }
  // Initializes the sequence object with a sequence in FASTA format read from
  // file stream.
  void read(FILE* in);
  // Prints the sequence in FASTA format to output stream.
  void write(FILE* fout, int width = 100) const;

  // Prints the Alignment in A2M format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const Sequence& seq) {
    const int kWidth = 100;
    out << '>' << seq.header_ << std::endl;
    for (int i = 0; i < seq.length(); ++i) {
      out << seq.chr(i);
      if ((i+1) % kWidth == 0) out << std::endl;
    }
    if (seq.length() % kWidth != 0) out << std::endl;
    return out;
  }

 private:
  // Buffer size for reading
  static const int kBufferSize = 16 * KB;

  // Convert the sequence in character representation to integer representation.
  void init(std::string header, std::string sequence);

  // The header without leading '>'.
  std::string header_;
  // The sequence itself in integer representation.
  std::valarray<char> seq_;
};  // Sequence

}  // namespace cs

#endif  // SRC_SEQUENCE_H_
