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
    double bit_score;
    // E-value of HSP
    double evalue;
    // Start of HSP in query
    int query_start;
    // End of HSP in query
    int query_end;
    // Query sequence segment
    std::vector<char> query;
    // Start of HSP in subject
    int subject_start;
    // End of HSP in subject
    int subject_end;
    // Query sequence segment
    std::vector<char> subject;
  };

  // Simple struct for data associated with each hit to a database sequence.
  struct Hit {
    // The ordinal id of the subject sequence this HSP list
    int oid;
    // Name of the subject sequence.
    std::string definition;
    // Smallest e-value for HSPs of this hit.
    double evalue;
    // Highest bit score for HSPs of this hit.
    double bit_score;
    // List of HSP's for one database sequence.
    std::vector<HSP> hsps;
  };

  typedef typename std::vector<Hit>::iterator HitIter;
  typedef typename std::vector<Hit>::const_iterator ConstHitIter;

  // Constructs results object from BLAST output in -m 0 format read from
  // file stream.
  explicit BlastResults(FILE* fin);

  ~BlastResults() {}

  // Returns a const iterator to the first hit.
  ConstHitIter begin() const { return hits_.begin(); }
  // Returns a const iterator just past the end of the hits vector.
  ConstHitIter end() const { return hits_.end(); }
  // Returns an iterator to the first integer element of the sequence.
  HitIter begin() { return &seq_[0]; }
  // Returns an iterator just past the end of the sequence.
  HitIter end() { return begin() + length(); }

  // Prints the results for logging
  friend std::ostream& operator<< (std::ostream& out, const Sequence& seq);

 private:
  // Initializes the sequence object with a sequence in FASTA format read from
  // file stream.
  void read(FILE* in);

};  // class BlastResults

}  // namespace cs

#endif  // SRC_BLAST_RESULTS_H_
