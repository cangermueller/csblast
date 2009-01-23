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
#include <valarray>
#include <vector>

#include "shared_ptr.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "matrix.h"
#include "exception.h"

namespace cs
{

class Alignment
{
  public:
    enum InputFormat {
        FASTA_INPUT    = 0,
        CLUSTAL_INPUT  = 1,
        A2M_INPUT      = 2,
        A3M_INPUT      = 2,
        PSI_INPUT      = 5,
        MAF_INPUT      = 6,
        BLAST_M6_INPUT = 7
    };
    enum OutputFormat {
        FASTA_OUTPUT   = 8,
        CLUSTAL_OUTPUT = 9,
        A2M_OUTPUT     = 10,
        A3M_OUTPUT     = 11,
        PSI_OUTPUT     = 12
    };

    typedef matrix<char>::row_type col_type;
    typedef matrix<char>::const_row_type const_col_type;

    // Constructs alignment multi FASTA formatted alignment read from input stream.
    Alignment(std::istream& in, InputFormat format, const SequenceAlphabet* alphabet);

    ~Alignment() {}

    // Access methods to get the integer representation of character in column i of sequence k.
    col_type operator[](int i) { return seqs_[match_cols_[i]]; }
    const_col_type operator[](int i) const { return seqs_[match_cols_[i]]; }
    // Returns the character in column i of sequence k.
    char chr(int k, int i) const { return alphabet_->itoc(seqs_[i][k]); }
    // Returns the number of sequences in the alignment.
    int nseqs() const { return seqs_.ncols(); }
    // Returns the total number of alignment columns.
    int ncols() const { return seqs_.nrows(); }
    // Returns the number of match columns.
    int nmatch() const { return match_cols_.size(); }
    // Returns the number of insert columns.
    int ninsert() const { return ncols() - nmatch(); }
    // Returns the header of sequence k.
    std::string header(int k) const { return headers_[k]; }
    // Sets the header of sequence k.
    void set_header(int k, const std::string& header) { headers_[k] = header; }
    // Remove all columns with a gap in the first sequence.
    void remove_columns_with_gap_in_first();
    // Remove all columns with more than X% gaps.
    void remove_columns_by_gap_rule(int gap_threshold = 50);
    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void read(std::istream& in, InputFormat format);
    // Writes the alignment in given format to ouput stream.
    void write(std::ostream& out, OutputFormat format, int width = 100) const;
    // Returns the underlying sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

  private:
    // Disallow copy and assign
    Alignment(const Alignment&);
    void operator=(const Alignment&);

    // Initializes alignment with given headers and sequences.
    void init(std::vector<std::string> headers, std::vector<std::string> seqs, bool case_insens = true);
    // Resize the sequence matrix and header vector to given dimensions. Attention: old data is lost!
    void resize(int nseqs, int ncols);
    // Reads an alignment in FASTA format read from given stream.
    void read_fasta(std::istream& in);
    // Writes the alignment in FASTA format to output stream.
    void write_fasta(std::ostream& out, int width = 100) const;

    // Row major matrix with sequences in integer representation.
    matrix<char> seqs_;
    // Array with indices of all columns [0,1,2,...,ncols-1].
    std::valarray<int> cols_;
    // Array with indices of match alignment columns.
    std::valarray<int> match_cols_;
    // Array mask indicating match and insert columns with true and false respectively.
    std::valarray<bool> mask_;
    // Headers of sequences in the alignment.
    std::vector<std::string> headers_;
    // Alphabet of sequences in the alignment.
    const SequenceAlphabet* alphabet_;
};  // Alignment



// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
// The returned pair vector specifies the sequence weights for each of the alignment sequences
// and the overall number of effective sequences in the alignment.
std::pair<std::vector<float>, float> global_weights_and_diversity(const Alignment& alignment);

// Calculates position-dependent sequence weights and number of effective sequences on subalignments.
// The return value is a pair consisting of a weights matrix (element (i,k) denotes the weight of
// sequence k in column i) and a vector with the number of effective sequences for alignment column.
std::pair< matrix<float>, std::vector<float> > position_specific_weights_and_diversity(const Alignment& alignment);

}//cs

#endif
