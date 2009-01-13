#ifndef CS_SEQUENCE_ALIGNMENT_H
#define CS_SEQUENCE_ALIGNMENT_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for an alignment of sequences over an alphabet.

#include <cctype>

#include <iostream>
#include <vector>

#include "smart_ptr.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "matrix.h"
#include "my_exception.h"

namespace cs
{

class Alignment
{
  public:
    friend std::istream& operator>> (std::istream& in, Alignment& alignment);

    // Constructs an alignment from alignment in multi FASTA format read from input stream.
    Alignment(std::istream& in, const SequenceAlphabet* alphabet);

    virtual ~Alignment() {}

    // Access methods to get the integer representation of character in column j of sequence i.
    char&       operator() (int i, int j) { return sequences_[i + j*nseqs_]; }
    const char& operator() (int i, int j) const { return sequences_[i + j*nseqs_]; }
    // Returns the character in column j of sequence i.
    char chr(int i, int j) const { return gap(i,j) ? kGap : alphabet_->itoc((*this)(i,j)); }
    // Returns true if the character at position (i,j) is a gap.
    bool gap(int i, int j) const { return (*this)(i,j) == gap(); }
    // Returns true if the character at position (i,j) is an endgap.
    bool endgap(int i, int j) const { return (*this)(i,j) == endgap(); }
    // Returns integer representation of gap.
    int gap() const { return alphabet_->size(); }
    // Returns integer representation of endgap.
    int endgap() const { return gap()+1; }
    // Returns the number of sequences in the alignment.
    int nseqs() const { return nseqs_; }
    // Returns the number of alignment columns.
    int ncols() const { return ncols_; }
    // Returns the header of sequence i.
    std::string header(int i) const { return headers_[i]; }
    // Sets the header of sequence i.
    void set_header(int i, const std::string& header) { headers_[i] = header; }
    // Returns the underlying sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

  private:
    // Gap character.
    static const char kGap = '-';

    // Disallow copy and assign
    Alignment(const Alignment&);
    void operator=(const Alignment&);

    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void init(std::istream& in);
    // Replaces gap with endgap for all gaps at either end of a sequence.
    void set_endgaps();
    // Resize the sequence matrix and header vector to given dimensions. Attention: old data is lost!
    void resize(int ncseqs, int ncols);

    // Number seqeuences in the alignment.
    int nseqs_;
    // Number alignment columns.
    int ncols_;
    // Row major matrix with sequences in integer representation.
    std::vector<char> sequences_;
    // Headers of sequences in the alignment.
    std::vector< std::string > headers_;
    // Alphabet of sequences in the alignment.
    const SequenceAlphabet* alphabet_;
};//Alignment



// Initializes a sequence alignment object from FASTA formatted alignment in input stream.
std::istream& operator>> (std::istream& in, Alignment& alignment);

// Prints the alignment in multi FASTA format to output stream.
std::ostream& operator<< (std::ostream& out, const Alignment& alignment);

// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
// The returned vector specifies the sequence weight for each of the alignment sequences.
std::vector<float> global_weights(const Alignment& alignment);

// Calculates position-dependent sequence weights and number of effective sequences on subalignments.
// The return value is a pair consisting of a weights matrix (element (i,j) denotes the weight of
// sequence j in column i) and a vector with the numbers of effective sequences in each alignment column.
std::pair< Matrix<float>, std::vector<float> > position_dependent_weights_and_neff(const Alignment& alignment);

}//cs

#endif
