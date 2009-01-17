#ifndef CS_ALIGNMENT_H
#define CS_ALIGNMENT_H
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
    friend std::ostream& operator<< (std::ostream& out, const Alignment& alignment);

    // Constructs alignment multi FASTA formatted alignment read from input stream.
    Alignment(std::istream& in, const SequenceAlphabet* alphabet);

    virtual ~Alignment() {}

    // Access methods to get the integer representation of character in column j of sequence i.
    char&       operator() (int i, int j) { return sequences_[i + j*nseqs_]; }
    const char& operator() (int i, int j) const { return sequences_[i + j*nseqs_]; }
    // Returns the character in column j of sequence i.
    char chr(int i, int j) const { return alphabet_->itoc((*this)(i,j)); }
     // Returns true if the character at position (i,j) is a real symbol (letter < ANY)
    bool less_any(int i, int j) const { return alphabet_->less_any((*this)(i,j)); }
    // Returns true if the character at position (i,j) is GAP.
    bool gap(int i, int j) const { return alphabet_->gap((*this)(i,j)); }
    // Returns true if the character at position (i,j) is ANY.
    bool any(int i, int j) const { return alphabet_->gap((*this)(i,j)); }
    // Returns true if the character at position (i,j) is ENDGAP.
    bool endgap(int i, int j) const { return alphabet_->endgap((*this)(i,j)); }
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

  protected:
    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    virtual void unserialize(std::istream& in);
    // Prints the alignment in multi FASTA format to ouput stream.
    virtual void serialize(std::ostream& out) const;

  private:
    // Disallow copy and assign
    Alignment(const Alignment&);
    void operator=(const Alignment&);

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
// The returned pair vector specifies the sequence weights for each of the alignment sequences
// and the overall number of effective sequences in the alignment.
std::pair<std::vector<float>, float> global_weights_and_diversity(const Alignment& alignment);

// Calculates position-dependent sequence weights and number of effective sequences on subalignments.
// The return value is a pair consisting of a weights matrix (element (i,k) denotes the weight of
// sequence k in column i) and a vector with the number of effective sequences for alignment column.
std::pair< Matrix<float>, std::vector<float> > position_dependent_weights_and_diversity(const Alignment& alignment);

}//cs

#endif
