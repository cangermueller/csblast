// Copyright 2009, Andreas Biegert

#ifndef SRC_BLAST_RESULTS_H_
#define SRC_BLAST_RESULTS_H_

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <string>
#include <vector>

#include "globals.h"
#include "utils.h"

namespace cs {

// A container class that class encapsulates all the search results of a
// BLAST search.
class BlastResults {
 public:
  // Simple struct for data associated with a HSP.
  struct HSP {
    // Number of identities in HSP
    int num_ident;
    // Score (in bits) of HSP
    float bit_score;
    // E-value of HSP
    double evalue;
    // Start of HSP in query
    int query_start;
    // End of HSP in query
    int query_end;
    // Query sequence segment
    std::vector<char> query_seg;
    // Start of HSP in subject
    int subject_start;
    // End of HSP in subject
    int subject_end;
    // Query sequence segment
    std::vector<char> subject_seg;
  };

  struct Hit {
    // The ordinal id of the subject sequence this HSP list
    int oid;
    // Smallest e-value for HSPs of this hit.
    double best_evalue;
    // List of HSP's for one database sequence.
    std::vector<HSP> hsps;
    // Name of the subject sequence.
    std::string subject_def;
  };

  typedef typename std::vector<Hit>::iterator HitIter;
  typedef typename std::vector<Hit>::const_iterator HitConstIter;

  // Constructs results object from BLAST output in -m 0 format read from
  // file stream.
  explicit BlastResults(FILE* fin);

  ~BlastResults() {}

  // Accessors for integer at position i of the sequence.
  Hit& operator[](int i) { return hits_[i]; }
  Hit char& operator[](int i) const { return hits_[i]; }
  Hit& hit(int i) { return hits_[i]; }
  const Hit& hit(int i) const { return hits_[i]; }

  // Returns a const iterator to the first hit.
  const_iterator begin() const { return hits_.begin(); }
  // Returns a const iterator just past the end of the hits vector.
  const_iterator end() const { return hits_.end(); }
  // Returns an iterator to the first integer element of the sequence.
  iterator begin() { return &seq_[0]; }
  // Returns an iterator just past the end of the sequence.
  iterator end() { return begin() + length(); }

  // Initializes the sequence object with a sequence in FASTA format read from
  // file stream.
  void read(FILE* in);
  // Prints the sequence in FASTA format to output stream.
  void write(FILE* fout, int width = 100) const;

  // Prints the Alignment in A2M format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const Sequence& seq) {
    const int kWidth = 100;
    out << '>' << seq.header_ << std::endl;
    for (int i = 0; i < seq.length(); ++i) {
      out << seq.chr(i);
      if ((i+1) % kWidth == 0) out << std::endl;
    }
    if (seq.length() % kWidth != 0) out << std::endl;
    return out;
  }

 private:
};  // class BlastResults

}  // namespace cs

#endif  // SRC_BLAST_RESULTS_H_
