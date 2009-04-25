// Copyright 2009, Andreas Biegert

#ifndef SRC_CSBLAST_ITERATION_H_
#define SRC_CSBLAST_ITERATION_H_

#include <string>
#include <vector>

#include "blast_hits.h"

namespace cs {

// Represents the iteration state in CS-BLAST
class CSBlastIteration {
 public:
  // List of sequence IDs
  typedef std::vector<std::string> SeqIds;

  // Number of iterations to perform. Use 0 to indicate that iterations must
  // take place until convergence.
  explicit CSBlastIteration(int num_iterations = 0);
  ~CSBlastIteration() {}

  // Returns the number of the current iteration
  int iteration() const { return iteration_; }
  // Returns true if search has converged (i.e.: no more new sequences have been
  // found since the last iteration)
  bool HasConverged() const;
  // Returns true if more iterations are still needed.
  bool HasMoreIterations() const;
  // Allows implicit conversion to a boolean value, returning true if there are
  // more iterations to perform or false if iterations are done.
  operator bool();
  // Advance the iterator by passing it the hits which passed the inclusion
  // criteria for the current iteration
  void Advance(const BlastHits& hits);

 private:
  // TODO: private data members
};  // class CSBlastIteration

}  // namespace cs

#endif  // SRC_CSBLAST_ITERATION_H_
