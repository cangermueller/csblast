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
    enum Format {
        FASTA    = 0,
        A2M      = 1,
        A3M      = 2,
        CLUSTAL  = 3,
        PSI      = 4
    };

    typedef matrix<char>::row_type col_type;
    typedef matrix<char>::const_row_type const_col_type;

    // Constructs alignment multi FASTA formatted alignment read from input stream.
    Alignment(std::istream& in, Format format, const SequenceAlphabet* alphabet);

    ~Alignment() {}

    // Access methods to get the integer representation of character in match column i of sequence k.
    col_type operator[](int i) { return seqs_[match_indexes[i]]; }
    const_col_type operator[](int i) const { return seqs_[match_indexes[i]]; }
    // Returns the integer in column i of sequence k.
    char seq(int k, int i) const { return seqs_[i][k]; }
    // Returns the character in column i of sequence k.
    char chr(int k, int i) const { return alphabet_->itoc(seqs_[i][k]); }
    // Returns the number of sequences in the alignment.
    int nseqs() const { return seqs_.ncols(); }
    // Returns the total number of alignment columns.
    int ncols() const { return seqs_.nrows(); }
    // Returns the number of match columns.
    int nmatch() const { return match_indexes.size(); }
    // Returns the number of insert columns.
    int ninsert() const { return ncols() - nmatch(); }
    // Returns the header of sequence k.
    std::string header(int k) const { return headers_[k]; }
    // Sets the header of sequence k.
    void set_header(int k, const std::string& header) { headers_[k] = header; }
    // Makes all columns with a residue in sequence k match columns.
    void assign_match_columns_by_sequence(int k = 0);
    // Makes all columns with less than X% gaps match columns.
    void assign_match_columns_by_gap_rule(int gap_threshold = 50);
    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void read(std::istream& in, Format format);
    // Writes the alignment in given format to ouput stream.
    void write(std::ostream& out, Format format, int width = 100) const;
    // Returns true if column i is a match column.
    bool match_column(int i) const { return match_column_[i]; }
    // Returns the underlying sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

  private:
    // Disallow copy and assign
    Alignment(const Alignment&);
    void operator=(const Alignment&);

    // Initializes alignment with given headers and sequences.
    void init(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);
    // Resize the sequence matrix and header vector to given dimensions. Attention: old data is lost!
    void resize(int nseqs, int ncols);
    // Fills match_indexes_ with the indexes of all match columns.
    void set_match_indexes();
    // Reads an alignment in FASTA format.
    void read_fasta(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Reads an alignment in A2M format from given stream.
    void read_a2m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Reads an alignment in A3M format from given stream.
    void read_a3m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
    void read_fasta_flavors(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Writes the alignment in FASTA, A2M, or A3M format to output stream.
    void write_fasta_flavors(std::ostream& out, Format format, int width = 100) const;

    // Row major matrix with sequences in integer representation.
    matrix<char> seqs_;
    // Array with indices of all columns [0,1,2,...,ncols-1].
    std::valarray<int> column_indexes_;
    // Array with indices of match columns.
    std::valarray<int> match_indexes;
    // Array mask indicating match and insert columns with true and false respectively.
    std::valarray<bool> match_column_;
    // Headers of sequences in the alignment.
    std::vector<std::string> headers_;
    // Alphabet of sequences in the alignment.
    const SequenceAlphabet* alphabet_;
};  // Alignment



// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
float global_weights_and_diversity(const Alignment& alignment, std::vector<float>& wg);

// Calculates position-dependent sequence weights and number of effective sequences on subalignments.
std::vector<float> position_specific_weights_and_diversity(const Alignment& alignment, matrix<float>& w);

}  // cs

#endif
