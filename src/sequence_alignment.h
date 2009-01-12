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

class SequenceAlignment
{
  public:
    friend std::istream& operator>> (std::istream& in, SequenceAlignment& alignment);

    SequenceAlignment(int nseqs,
                      int ncols,
                      const SequenceAlphabet* alphabet);
    SequenceAlignment(std::istream& in,
                      const SequenceAlphabet* alphabet);
    virtual ~SequenceAlignment();

    // Access methods to get the integer representation in column j of sequence i (starting at index 0)
    char&       operator() (int i, int j) { return sequences_[i + j*nseqs_]; }
    const char& operator() (int i, int j) const { return sequences_[i + j*nseqs_]; }
    // Returns the character in column j of sequence i (starting at index 0)
    char chr(int i, int j) const { return gap(i,j) ? kGap : alphabet_->itoc((*this)(i,j)); }
    // Returns true if the character at position (i,j) is a gap.
    bool gap(int i, int j) const { return (*this)(i,j) == gaptoi(); }
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
    // Gap character
    static const char kGap = '-';

    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void init(std::istream& in);
    // Returns integer representation of a gap.
    int gaptoi() const { return alphabet_->size(); }
    // Resize the sequence matrix and header vector to given dimensions. Attention: old data is lost!
    void resize(int ncseqs, int ncols);

    // Number seqeuences in the alignment
    int nseqs_;
    // Number alignment columns
    int ncols_;
    // Row major matrix with sequences in integer representation
    std::vector<char> sequences_;
    // Headers of sequences in the alignment
    std::vector< std::string > headers_;
    // Alphabet of sequences in the alignment
    const SequenceAlphabet* alphabet_;
};//SequenceAlignment



// Initializes a sequence alignment object from FASTA formatted alignment in input stream.
std::istream& operator>> (std::istream& in, SequenceAlignment& alignment);

// Prints the alignment in multi FASTA format to output stream.
std::ostream& operator<< (std::ostream& out, const SequenceAlignment& alignment);

// Calculates column specific sequence weights from subalignments within the global alignment.
Matrix<float> column_specific_sequence_weights(const SequenceAlignment& alignment);

// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
std::vector<float> global_sequence_weights(const SequenceAlignment& alignment);

}//cs

#endif
