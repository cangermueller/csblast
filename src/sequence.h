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
#include <vector>

#include "sequence_alphabet.h"
#include "smart_ptr.h"

namespace cs
{

class Sequence
{
  public:
    typedef std::vector<char>::const_iterator const_iterator;
    typedef std::vector<char>::iterator iterator;

    friend std::istream& operator>> (std::istream& in, Sequence& sequence);

    // Constructs sequence with specified length and alphabet.
    Sequence(int length, const SequenceAlphabet* alphabet);
    // Constructs sequence from serialized sequence in FASTA format read from input stream.
    Sequence(std::istream& in, const SequenceAlphabet* alphabet);
    // Constructs sequence with given header and sequence string of characters.
    Sequence(const std::string& header,
             const std::string& sequence,
             const SequenceAlphabet* alphabet);

    virtual ~Sequence() {}

    // Reads all available sequences from the input stream and returns them in a vector.
    static std::vector< SmartPtr<Sequence> > read(std::istream& in,
                                                  const SequenceAlphabet* alphabet);
    // Accessors for character at position i of the sequence in integer representation.
    char&       operator() (int i) { return sequence_[i]; }
    const char& operator() (int i) const { return sequence_[i]; }
    // Returns the character at position i of the sequence.
    char chr(int i) const { return alphabet_->itoc(sequence_[i]); }
    // Returns the sequence length.
    int length() const { return sequence_.size(); }
    // Returns the header string of the sequence.
    std::string header() const { return header_; }
    // Sets the header to given string.
    void set_header(const std::string& header) { header_ = header; }
    // Returns the sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }
    // Returns a const iterator to the first integer element of the sequence.
    const_iterator begin() const { return sequence_.begin(); }
    // Returns a const iterator just past the end of the sequence.
    const_iterator end() const { return sequence_.end(); }
    // Returns an iterator to the first integer element of the sequence.
    iterator begin() { return sequence_.begin(); }
    // Returns an iterator just past the end of the sequence.
    iterator end() { return sequence_.begin(); }

  private:
    // Disallow copy and assign
    Sequence(const Sequence&);
    void operator=(const Sequence&);

    // Initializes the sequence object with a sequence in FASTA format read from given stream.
    void init(std::istream& in);
    // Convert the sequence in character representation to integer representation.
    void check_and_convert();

    const SequenceAlphabet* alphabet_;
    std::string header_;
    std::vector<char> sequence_;
};//Sequence



// Initializes a sequence object from FASTA formatted sequence in input stream.
std::istream& operator>> (std::istream& in, Sequence& sequence);

// Prints the sequence in FASTA format to output stream.
std::ostream& operator<< (std::ostream& out, const Sequence& sequence);

}//cs

#endif
