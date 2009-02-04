#ifndef CS_SEQUENCE_H
#define CS_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class representing a sequence consisting of letters over a
// sequence alphabet.

#include <cctype>

#include <algorithm>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "exception.h"
#include "log.h"
#include "shared_ptr.h"

namespace cs
{

template<class Alphabet_T>
class Sequence
{
  public:
    typedef char* iterator;
    typedef const char* const_iterator;

    // Constructs sequence with specified length.
    Sequence(int length);
    // Constructs sequence from serialized sequence in FASTA format read from input stream.
    Sequence(std::istream& in);
    // Constructs sequence with given header and sequence string of characters.
    Sequence(const std::string& header,
             const std::string& sequence);

    ~Sequence() {}

    // Reads all available sequences from the input stream and returns them in a vector.
    static std::vector< shared_ptr<Sequence> > readall(std::istream& in);

    // Accessors for integer at position i of the sequence.
    char& operator[](int i) { return seq_[i]; }
    const char& operator[](int i) const { return seq_[i]; }
    // Returns the character at position i of the sequence.
    char chr(int i) const { return Alphabet_T::instance().itoc(seq_[i]); }
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
    // Initializes the sequence object with a sequence in FASTA format read from given stream.
    void read(std::istream& in);
    // Prints the sequence in FASTA format to output stream.
    void write(std::ostream& out, int width = 100) const;

    // Prints the Alignment in A2M format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const Sequence& sequence)
    {
        sequence.write(out);
        return out;
    }

  private:
    // Disallow copy and assign
    Sequence(const Sequence&);
    void operator=(const Sequence&);

    // Convert the sequence in character representation to integer representation.
    void init(std::string header, std::string sequence);

    // The header without leading '>'.
    std::string header_;
    // The sequence itself in integer representation.
    std::valarray<char> seq_;
};  // Sequence


template<class Alphabet_T>
Sequence<Alphabet_T>::Sequence(int length)
        : seq_(length)
{}

template<class Alphabet_T>
Sequence<Alphabet_T>::Sequence(std::istream& in)
{ read(in); }

template<class Alphabet_T>
Sequence<Alphabet_T>::Sequence(const std::string& header,
                   const std::string& sequence)
{ init(header, sequence); }

template<class Alphabet_T>
std::vector< shared_ptr< Sequence<Alphabet_T> > > Sequence<Alphabet_T>::readall(std::istream& in)
{
    std::vector< shared_ptr<Sequence> > sequences;
    while (in.good()) {
        shared_ptr<Sequence> p(new Sequence(in));
        sequences.push_back(p);
    }

    return sequences;
}

template<class Alphabet_T>
void Sequence<Alphabet_T>::init(std::string header, std::string sequence)
{
    // init header
    std::string().swap(header_);
    header_.append(header.begin() + (header[0]=='>' ? 1 : 0), header.end());

    //strip whitespace and newlines from sequence.
    sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace), sequence.end());
    //validate each character and convert to integer representation
    const Alphabet_T& alphabet = Alphabet_T::instance();
    const int seqlen = sequence.length();
    seq_.resize(seqlen);
    for (int i = 0; i < seqlen; ++i) {
        char c = sequence[i];
        if (alphabet.valid(c)) {
            seq_[i] = alphabet.ctoi(c);
        } else {
            throw Exception("Invalid character %c at position %i of sequence '%s'", c, i, header_.c_str());
        }
    }
}

template<class Alphabet_T>
void Sequence<Alphabet_T>::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading sequence from stream ...";
    std::string tmp;
    std::string header;
    std::string sequence;

    //read header
    if (getline(in, tmp)) {
        if (tmp.empty() ||  tmp[0] != '>')
            throw Exception("Bad format: first line of FASTA sequence does not start with '>' character!");
        header.append(tmp.begin() + 1, tmp.end());
    } else {
        throw Exception("Failed to read from FASTA formatted input stream!");
    }
    //read sequence
    while (in.peek() != '>' && getline(in, tmp))
        sequence.append(tmp.begin(), tmp.end());

    init(header, sequence);
    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void Sequence<Alphabet_T>::write(std::ostream& out, int width) const
{
    out << '>' << header_ << std::endl;
    for (int i = 0; i < length(); ++i) {
        out << chr(i);
        if ((i+1) % width == 0) out << std::endl;
    }
    if (length() % width != 0) out << std::endl;
}

}  // cs

#endif
