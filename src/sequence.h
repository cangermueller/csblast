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

    Sequence(int length,
             const SequenceAlphabet* alphabet);
    Sequence(const std::string& header,
             const std::string& sequence,
             const SequenceAlphabet* alphabet);
    Sequence(std::istream& in,
             const SequenceAlphabet* alphabet);
    virtual ~Sequence();

    // Reads all available sequences from the input stream and returns them in a vector.
    static std::vector< SmartPtr<Sequence> > read(std::istream& in,
                                                  const SequenceAlphabet* alphabet);
    char&       operator() (int i) { return sequence_[i]; }
    const char& operator() (int i) const { return sequence_[i]; }
    char chr(int i) const { return alphabet_->itoc(sequence_[i]); }
    int length() const { return sequence_.size(); }
    std::string header() const { return header_; }
    void set_header(const std::string& header) { header_ = header; }
    const SequenceAlphabet* alphabet() const { return alphabet_; }
    const_iterator begin() const { return sequence_.begin(); }
    const_iterator end() const { return sequence_.end(); }
    iterator begin() { return sequence_.begin(); }
    iterator end() { return sequence_.begin(); }

  private:
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
