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
#include "blast_hits.h"
#include "matrix.h"
#include "sequence.h"
#include "shared_ptr.h"

namespace cs {

// Forward declarations
template<class Alphabet>
class Alignment;

// Convince the compiler that operator<< is a template friend.
template<class Alphabet>
std::ostream& operator<< (std::ostream& out, const Alignment<Alphabet>& ali);

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

  // Constructs alignment from multi FASTA formatted alignment read from input
  // stream.
  Alignment(FILE* fin, Format format);
  // Constructs an alignment from a single sequence.
  Alignment(const Sequence<Alphabet>& seq);
  // Constructs a query anchored alignment from BLAST hits. If flag best is set
  // to true only the best HSP of each hit is included in the alignment.
  Alignment(const BlastHits& hits, bool best = false);

  ~Alignment() {}

  // Reads all available alignments from the input stream and returns them in a
  // vector.
  static void ReadAll(FILE* fin,
                      Format format,
                      std::vector< shared_ptr<Alignment> >* v);

  // Accessors for integer representation of character in MATCH column i of
  // sequence k.
  col_type operator[](int i) { return seqs_[match_indices[i]]; }
  const_col_type operator[](int i) const { return seqs_[match_indices[i]]; }
  char& at(int i, int k) { return seqs_[match_indices[i]][k]; }
  const char& at(int i, int k) const { return seqs_[match_indices[i]][k]; }
  // Accessors for integer representation of the character at position i
  // (NOT match column i) of sequence k.
  char& seq(int k, int i) { return seqs_[i][k]; }
  const char& seq(int k, int i) const { return seqs_[i][k]; }
  // Returns the character at position i (NOT match column i) of sequence k.
  char chr(int k, int i) const {
    return Alphabet::instance().itoc(seqs_[i][k]);
  }
  // Returns the number of sequences in the alignment.
  int num_seqs() const { return seqs_.num_cols(); }
  // Returns the total number of alignment columns.
  int num_cols() const { return seqs_.num_rows(); }
  // Returns the number of match columns.
  int num_match_cols() const { return match_indices.size(); }
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
  void AssignMatchColumnsBySequence(int k = 0);
  // Makes all columns with less than X% gaps match columns.
  void AssignMatchColumnsByGapRule(int gap_threshold = 50);
  // Initializes object with an alignment in FASTA format read from given
  // stream.
  void Read(FILE* fin, Format format);
  // Writes the alignment in given format to ouput stream.
  void Write(FILE* fout, Format format, int width = 100) const;
  // Returns true if column i is a match column.
  bool IsMatchColumn(int i) const { return match_column_[i]; }
  // Removes all insert columns from the alignment.
  void RemoveInsertColumns();
  // Merges the provided alignment with this alignment considering only sequences
  // that are not already included in this alignment. Warning: Inserts in current
  // alignment are lost!
  void Merge(const Alignment<Alphabet>& ali);
 // Returns the sequence k as Sequence object.
  Sequence<Alphabet> GetSequence(int k) const;

  // Prints the Alignment in A2M format for debugging.
  friend std::ostream& operator<< <> (std::ostream& out,
                                      const Alignment<Alphabet>& ali);

 private:
  // Buffer size for reading
  static const int kBufferSize = 32 * KB;

  // Initializes alignment with given headers and sequences.
  void Init(const std::vector<std::string>& headers,
            const std::vector<std::string>& seqs);
  // Resize the sequence matrix and header vector to given dimensions.
  void Resize(int num_seqs, int num_cols);
  // Fills match_indices_ with the indices of all match columns.
  void SetMatchIndices();
  // Reads an alignment in FASTA format.
  void ReadFasta(FILE* fin,
                 std::vector<std::string>* headers,
                 std::vector<std::string>* seqs);
  // Reads an alignment in A2M format from given stream.
  void ReadA2M(FILE* fin,
               std::vector<std::string>* headers,
               std::vector<std::string>* seqs);
  // Reads an alignment in A3M format from given stream.
  void ReadA3M(FILE* fin,
               std::vector<std::string>* headers,
               std::vector<std::string>* sequences);
  // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
  void ReadFastaFlavors(FILE* fin,
                        std::vector<std::string>* headers,
                        std::vector<std::string>* seqs);
  // Reads an alignment in PSI format.
  void ReadPsi(FILE* fin,
               std::vector<std::string>* headers,
               std::vector<std::string>* seqs);
  // Writes the alignment in FASTA, A2M, or A3M format to output stream.
  void WriteFastaFlavors(FILE* fout, Format format, int width = 100) const;
  // Writes the alignment in CLUSTAL or PSI format to output stream.
  void WriteClustalFlavors(FILE* fout, Format format, int width = 100) const;

  // Row major matrix with sequences in integer representation.
  matrix<char> seqs_;
  // Array with indices of all columns [0,1,2,...,num_cols-1].
  std::valarray<int> column_indices_;
  // Array with indices of match columns.
  std::valarray<int> match_indices;
  // Array mask indicating match and insert columns.
  std::valarray<bool> match_column_;
  // Headers of sequences in the alignment.
  std::vector<std::string> headers_;
};  // Alignment



// Calculates global sequence weights by maximum entropy weighting
// (Henikoff&Henikoff '94).
template<class Alphabet>
float GlobalWeightsAndDiversity(const Alignment<Alphabet>& alignment,
                                   std::vector<float>& wg);

// Calculates position-dependent sequence weights and number of effective
// sequences on subalignments.
template<class Alphabet>
std::vector<float> PositionSpecificWeightsAndDiversity(
    const Alignment<Alphabet>& alignment, matrix<float>& w);

// Returns the alignment format corresponding to provided filename extension
template<class Alphabet>
typename Alignment<Alphabet>::Format AlignmentFormatFromString(
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
