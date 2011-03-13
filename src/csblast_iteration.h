// Copyright 2009, Andreas Biegert

#ifndef CS_CSBLAST_ITERATION_H_
#define CS_CSBLAST_ITERATION_H_

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
  int IterationNumber() const { return iterations_done_ + 1; }

  // Allows implicit conversion to a boolean value, returning true if there are
  // more iterations to perform or false if iterations are done.
  operator bool() const { return (HasMoreIterations() && !HasConverged()); }

  // Returns true if more iterations are still needed.
  bool HasMoreIterations() const {
    return (iterations_todo_ == 0 || iterations_done_ < iterations_todo_);
  }

  // Returns true if search has converged (i.e.: no more new sequences have been
  // found since the last iteration)
  bool HasConverged() const;

  // Advance the iterator by passing it the hits which passed the inclusion
  // criteria for the current iteration
  void Advance(const BlastHits& hits);

 private:
  // Number of iterations to perform
  int iterations_todo_;
  // Number of iterations already performed
  int iterations_done_;
  // Identifiers for sequences found in the previous iteration
  SeqIds previous_data_;
  // Identifiers for sequences found in the current iteration
  SeqIds current_data_;
};  // class CSBlastIteration

}  // namespace cs

#endif  // CS_CSBLAST_ITERATION_H_
