#ifndef CS_SEQUENCE_H
#define CS_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class representing a sequence consisting of letters over a
// sequence alphabet.

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <streambuf>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "sequence_alphabet.h"
#include "my_exception.h"

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
    Sequence(std::string::const_iterator start,
             std::string::const_iterator end,
             const SequenceAlphabet* alphabet);
    virtual ~Sequence();

    static std::vector<Sequence*> read(std::istream& in,
                                       const SequenceAlphabet* alphabet);

    char&       operator() (int i);
    const char& operator() (int i) const;
    char chr(int i) const;
    int length() const;
    const std::string& header() const;
    std::string& header();
    const SequenceAlphabet& alphabet() const;
    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

private:
    // Initializes the sequence object with a header and a sequence of characters.
    void init(const std::string& header, const std::string& sequence);

    const SequenceAlphabet* alphabet_;
    std::string header_;
    std::vector<char> sequence_;
};//Sequence


// Initializes a sequence object from FASTA formatted sequence in input stream.
std::istream& operator>> (std::istream& in, Sequence& sequence);

// Prints the sequence in FASTA format to output stream.
std::ostream& operator<< (std::ostream& out, const Sequence& sequence);

// Parse header and sequence from data starting from given index. After the parse index points to the start of
// the next sequence or is zero if there is none.
std::pair<std::string, std::string > parse_fasta_sequence(const std::string& data, size_t& index);


inline char& Sequence::operator() (int i)
{ return sequence_[i]; }

inline const char& Sequence::operator() (int i) const
{ return sequence_[i]; }

inline char Sequence::chr(int i) const
{ return alphabet_->itoc(sequence_[i]); }

inline int Sequence::length() const
{ return sequence_.size(); }

inline std::string& Sequence::header()
{ return header_; }

inline const std::string& Sequence::header() const
{ return header_; }

inline const SequenceAlphabet& Sequence::alphabet() const
{ return *alphabet_; }

inline Sequence::const_iterator Sequence::begin() const
{ return sequence_.begin(); }

inline Sequence::const_iterator Sequence::end() const
{ return sequence_.end(); }

inline Sequence::iterator Sequence::begin()
{ return sequence_.begin(); }

inline Sequence::iterator Sequence::end()
{ return sequence_.end(); }

}//cs

#endif
