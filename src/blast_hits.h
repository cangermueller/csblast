// Copyright 2009, Andreas Biegert

#ifndef SRC_BLAST_HITS_H_
#define SRC_BLAST_HITS_H_

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <string>
#include <vector>

#include "globals.h"
#include "utils-inl.h"

namespace cs {

// A container class that class encapsulates all the search results of a
// BLAST search.
class BlastHits {
 public:
  // Simple struct for data associated with a HSP.
  struct Hsp {
    Hsp() : bit_score(0.0),
            evalue(0.0),
            query_start(0),
            query_end(0),
            subject_start(0),
            subject_end(0) {}

    // Score (in bits) of HSP
    double bit_score;
    // E-value of HSP
    double evalue;
    // Start of HSP in query
    int query_start;
    // End of HSP in query
    int query_end;
    // Query sequence segment
    std::vector<char> query_seq;
    // Start of HSP in subject
    int subject_start;
    // End of HSP in subject
    int subject_end;
    // Query sequence segment
    std::vector<char> subject_seq;
    // Length of the HSP alignment
    int length;
  };

  // Simple struct for data associated with each hit to a database sequence.
  struct Hit {
    Hit() : oid(0), evalue(0.0), bit_score(0.0) {}

    // The ordinal id of the subject sequence this HSP list
    int oid;
    // Name of the subject sequence.
    std::string definition;
    // Smallest e-value for HSPs of this hit.
    double evalue;
    // Highest bit score for HSPs of this hit.
    double bit_score;
    // List of HSP's for one database sequence.
    std::vector<Hsp> hsps;
  };

  typedef std::vector<Hit>::iterator HitIter;
  typedef std::vector<Hit>::const_iterator ConstHitIter;
  typedef std::vector<Hsp>::iterator HspIter;
  typedef std::vector<Hsp>::const_iterator ConstHspIter;

  // Constructs an empty hits object that can be filled by calling Read
  BlastHits();
  // Constructs an hit object from BLAST output in -m 0 format read from
  // file stream.
  explicit BlastHits(FILE* fin);

  ~BlastHits() {}

  // Accessors for integer at position i of the sequence.
  Hit& operator[](int i) { return hits_[i]; }
  const Hit& operator[](int i) const { return hits_[i]; }
  Hit& hit(int i) { return hits_[i]; }
  const Hit& hit(int i) const { return hits_[i]; }
  // Returns a const iterator to the first hit.
  ConstHitIter begin() const { return hits_.begin(); }
  // Returns a const iterator just past the end of the hits vector.
  ConstHitIter end() const { return hits_.end(); }
  // Returns an iterator to the first integer element of the sequence.
  HitIter begin() { return hits_.begin(); }
  // Returns an iterator just past the end of the sequence.
  HitIter end() { return hits_.end(); }
  // Returns number of hits.
  int num_hits() const { return hits_.size(); }
  // Returns number of hits.
  int size() const { return hits_.size(); }
  // Returns length of query sequence.
  int query_length() const { return query_length_; }
  // Returns true if hit list is empty.
  bool empty() const { return hits_.empty(); }
  // Filters hits by e-value threshold.
  void Filter(double evalue_threshold);
  // Fills the hits object with with hits parsed from BLAST output.
  void Read(FILE* fin);

  // Prints the results for logging
  friend std::ostream& operator<< (std::ostream& out, const BlastHits& res);

 private:
  // List of hits in the BLAST results
  std::vector<Hit> hits_;
  // Length of query sequence
  int query_length_;
};  // class BlastHits

}  // namespace cs

#endif  // SRC_BLAST_HITS_H_
