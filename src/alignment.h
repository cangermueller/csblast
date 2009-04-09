// Copyright 2009, Andreas Biegert

#ifndef SRC_ALIGNMENT_H_
#define SRC_ALIGNMENT_H_

#include <cctype>
#include <cstdio>

#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "globals.h"
#include "log.h"
#include "matrix.h"
#include "exception.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs {

// A container class for multiple sequence alignments.
template<class Alphabet>
class Alignment {
 public:
  // Supported alignment formats for in- and output.
  enum Format {
    FASTA    = 0,
    A2M      = 1,
    A3M      = 2,
    CLUSTAL  = 3,
    PSI      = 4
  };

  // Public typedefs
  typedef matrix<char>::row_type col_type;
  typedef matrix<char>::const_row_type const_col_type;
  typedef matrix<char>::col_type row_type;
  typedef matrix<char>::const_col_type const_row_type;
  typedef matrix<char>::iterator iterator;
  typedef matrix<char>::const_iterator const_iterator;

  // Constructs alignment multi FASTA formatted alignment read from input
  // stream.
  Alignment(std::istream& in, Format format);
  // Constructs alignment multi FASTA formatted alignment read from input
  // stream.
  Alignment(FILE* fin, Format format);

  ~Alignment() {}

  // Reads all available alignments from the input stream and returns them in a
  // vector.
  static void readall(std::istream& in,
                      Format format,
                      std::vector< shared_ptr<Alignment> >& v);
  // Reads all available alignments from the input stream and returns them in a
  // vector.
  static void readall(FILE* fin,
                      Format format,
                      std::vector< shared_ptr<Alignment> >* v);

  // Accessors for integer representation of character in match column i of
  // sequence k.
  col_type operator[](int i) { return seqs_[match_indexes[i]]; }
  const_col_type operator[](int i) const { return seqs_[match_indexes[i]]; }
  // Returns the integer in column i of sequence k.
  char seq(int k, int i) const { return seqs_[i][k]; }
  // Returns the character in column i of sequence k.
  char chr(int k, int i) const {
    return Alphabet::instance().itoc(seqs_[i][k]);
  }
  // Returns the number of sequences in the alignment.
  int num_seqs() const { return seqs_.num_cols(); }
  // Returns the total number of alignment columns.
  int num_cols() const { return seqs_.num_rows(); }
  // Returns the number of match columns.
  int num_match_cols() const { return match_indexes.size(); }
  // Returns the number of insert columns.
  int num_insert_cols() const { return num_cols() - num_match_cols(); }
  // Returns the total number of characters in the alignment (incl. inserts).
  int size() const { return seqs_.size(); }

  row_type seq_begin(int k) { return seqs_.col_begin(k); }
  // Returns an iterator just past the end of alignment sequence k.
  row_type seq_end(int k) { return seqs_.col_end(k); }
  // // Returns a const iterator to the first element in alignment sequence k.
  const_row_type seq_begin(int k) const { return seqs_.col_begin(k); }
  // // Returns a const iterator just past the end of alignment sequence k.
  const_row_type seq_end(int k) const { return seqs_.col_end(k); }

  // Returns an iterator to the first element in the alignment matrix
  // (incl. inserts).
  iterator begin() { return seqs_.begin(); }
  // Returns an iterator just past the end of the alignment matrix
  // (incl. inserts).
  iterator end() { return seqs_.end(); }
  // Returns a const iterator to the first element in the alignment matrix
  // (incl. inserts).
  const_iterator begin() const { return seqs_.begin(); }
  // Returns a const iterator just past the end of the alignment matrix
  // (incl. inserts).
  const_iterator end() const { return seqs_.end(); }

  // Returns the header of sequence k.
  std::string header(int k) const { return headers_[k]; }
  // Sets the header of sequence k.
  void set_header(int k, const std::string& header) { headers_[k] = header; }
  // Makes all columns with a residue in sequence k match columns.
  void assign_match_columns_by_sequence(int k = 0);
  // Makes all columns with less than X% gaps match columns.
  void assign_match_columns_by_gap_rule(int gap_threshold = 50);
  // Initializes object with an alignment in FASTA format read from given
  // stream.
  void read(std::istream& in, Format format);
  // Initializes object with an alignment in FASTA format read from given
  // stream.
  void read(FILE* fin, Format format);
  // Writes the alignment in given format to ouput stream.
  void write(std::ostream& out, Format format, int width = 100) const;
  // Returns true if column i is a match column.
  bool match_column(int i) const { return match_column_[i]; }
  // Removes all insert columns from the alignment.
  void remove_insert_columns();

  // Prints the Alignment in A2M format for debugging.
  friend std::ostream& operator<< (std::ostream& out,
                                   const Alignment& alignment) {
    alignment.write(out, Alignment::CLUSTAL);
    return out;
  }

 private:
  // Buffer size for reading
  static const int kBufferSize = 32 * KB;

  // Disallow copy and assign
  Alignment(const Alignment&);
  void operator=(const Alignment&);

  // Initializes alignment with given headers and sequences.
  void init(const std::vector<std::string>& headers,
            const std::vector<std::string>& seqs);
  // Resize the sequence matrix and header vector to given dimensions.
  void resize(int num_seqs, int num_cols);
  // Fills match_indexes_ with the indexes of all match columns.
  void set_match_indexes();
  // Reads an alignment in FASTA format.
  void read_fasta(std::istream& in,
                  std::vector<std::string>& headers,
                  std::vector<std::string>& seqs);
  // Reads an alignment in A2M format from given stream.
  void read_a2m(std::istream& in,
                std::vector<std::string>& headers,
                std::vector<std::string>& seqs);
  // Reads an alignment in A3M format from given stream.
  void read_a3m(std::istream& in,
                std::vector<std::string>& headers,
                std::vector<std::string>& seqs);
  // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
  void read_fasta_flavors(std::istream& in,
                          std::vector<std::string>& headers,
                          std::vector<std::string>& seqs);
  // Reads an alignment in FASTA format.
  void read_fasta(FILE* fin,
                  std::vector<std::string>* headers,
                  std::vector<std::string>* seqs);
  // Reads an alignment in A2M format from given stream.
  void read_a2m(FILE* fin,
                std::vector<std::string>* headers,
                std::vector<std::string>* seqs);
  // Reads an alignment in A3M format from given stream.
  void read_a3m(FILE* fin,
                std::vector<std::string>* headers,
                std::vector<std::string>* sequences);
  // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
  void read_fasta_flavors(FILE* fin,
                          std::vector<std::string>* headers,
                          std::vector<std::string>* seqs);
  // Writes the alignment in FASTA, A2M, or A3M format to output stream.
  void write_fasta_flavors(std::ostream& out,
                           Format format,
                           int width = 100) const;
  // Writes the alignment in CLUSTAL or PSI format to output stream.
  void write_clustal_flavors(std::ostream& out,
                             Format format,
                             int width = 100) const;

  // Row major matrix with sequences in integer representation.
  matrix<char> seqs_;
  // Array with indices of all columns [0,1,2,...,num_cols-1].
  std::valarray<int> column_indexes_;
  // Array with indices of match columns.
  std::valarray<int> match_indexes;
  // Array mask indicating match and insert columns.
  std::valarray<bool> match_column_;
  // Headers of sequences in the alignment.
  std::vector<std::string> headers_;
};  // Alignment



// Calculates global sequence weights by maximum entropy weighting
// (Henikoff&Henikoff '94).
template<class Alphabet>
float global_weights_and_diversity(const Alignment<Alphabet>& alignment,
                                   std::vector<float>& wg);

// Calculates position-dependent sequence weights and number of effective
// sequences on subalignments.
template<class Alphabet>
std::vector<float> position_specific_weights_and_diversity(
    const Alignment<Alphabet>& alignment, matrix<float>& w);

// Returns the alignment format corresponding to provided filename extension
template<class Alphabet>
typename Alignment<Alphabet>::Format alignment_format_from_string(
    const std::string& s);

// Converts a character to uppercase and '.' to '-'.
inline char to_match_chr(char c) {
  return isalpha(c) ? toupper(c) : (c == '.' ? '-' : c);
}

// Converts a character to lowercase and '-' to '.'.
inline char to_insert_chr(char c) {
  return isalpha(c) ? tolower(c) : (c == '-' ? '.' : c);
}

// Predicate indicating if character belongs to match column.
inline bool match_chr(char c) {
  return (isalpha(c) && isupper(c)) || c == '-';
}

// Predicate indicating if character belongs to insert column.
inline char insert_chr(char c) {
  return (isalpha(c) && islower(c)) || c == '.';
}

}  // namespace cs

#endif  // SRC_ALIGNMENT_H_
