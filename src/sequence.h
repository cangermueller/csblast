#ifndef CS_SEQUENCE_H
#define CS_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class representing a sequence consisting of letters over a
// sequence alphabet.

#include <string>
#include <valarray>
#include <vector>

#include "sequence_alphabet.h"
#include "shared_ptr.h"

namespace cs
{

class Sequence
{
  public:
    typedef char* iterator;
    typedef const char* const_iterator;

    // Constructs sequence with specified length and alphabet.
    Sequence(int length, const SequenceAlphabet* alphabet);
    // Constructs sequence from serialized sequence in FASTA format read from input stream.
    Sequence(std::istream& in, const SequenceAlphabet* alphabet);
    // Constructs sequence with given header and sequence string of characters.
    Sequence(const std::string& header,
             const std::string& sequence,
             const SequenceAlphabet* alphabet);

    ~Sequence() {}

    // Reads all available sequences from the input stream and returns them in a vector.
    static std::vector< shared_ptr<Sequence> > readall(std::istream& in,
                                                       const SequenceAlphabet* alphabet);

    // Accessors for integer at position i of the sequence.
    char& operator[](int i) { return seq_[i]; }
    const char& operator[](int i) const { return seq_[i]; }
    // Returns the character at position i of the sequence.
    char chr(int i) const { return alphabet_->itoc(seq_[i]); }
    // Returns the sequence length.
    int length() const { return seq_.size(); }
    // Returns the header string of the sequence.
    std::string header() const { return header_; }
    // Sets the header to given string.
    void set_header(const std::string& header) { header_ = header; }
    // Returns the sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }
    // Returns a const iterator to the first integer element of the sequence.
    const_iterator begin() const { return &seq_[0]; }
    // Returns a const iterator just past the end of the sequence.
    const_iterator end() const { return begin() + length(); }
    // Returns an iterator to the first integer element of the sequence.
    iterator begin() { return &seq_[0]; }
    // Returns an iterator just past the end of the sequence.
    iterator end() { return begin() + length(); }
    // Initializes the sequence object with a sequence in FASTA format read from given stream.
    void read(std::istream& in);
    // Prints the sequence in FASTA format to output stream.
    void write(std::ostream& out, int width = 100) const;

  private:
    // Disallow copy and assign
    Sequence(const Sequence&);
    void operator=(const Sequence&);

    // Convert the sequence in character representation to integer representation.
    void init(std::string header, std::string sequence);

    // The underlying alphabet of the sequence.
    const SequenceAlphabet* alphabet_;
    // The header without leading '>'.
    std::string header_;
    // The sequence itself in integer representation.
    std::valarray<char> seq_;
};  // Sequence

}//cs

#endif
