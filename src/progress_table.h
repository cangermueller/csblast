// Copyright 2009, Andreas Biegert

#ifndef CS_PROGRESS_TABLE_H_
#define CS_PROGRESS_TABLE_H_

#include <stdlib.h>
#include <stdio.h>

namespace cs {

// Abstract base class for progress table to be subclassed for HMM training
// and EM clustering.
class ProgressTable {
 public:
  // To be used by derived classes.
  ProgressTable(FILE* fout = stdout, int width = 30);

  virtual ~ProgressTable() {}

  // Prints header information.
  virtual void PrintHeader() = 0;
  // Starts a new row and prints statistics to outstream.
  virtual void PrintRowBegin() = 0;
  // Ends the current row and prints likelihood.
  virtual void PrintRowEnd() = 0;
  // Sets total work per scan.
  void set_total_work(long work) { work_total_ = work; }
  // Resets the progress bar to zero.
  void Reset();
  // Advances the progress bar proportional to the amount of work done.
  void print_progress(long work);
  // Returns output stream of the progress table
  FILE* stream() const { return fout_; }

 protected:

  // Output stream.
  FILE* fout_;
  // Width of prograess bar.
  const int width_;
  // Currnt width of bar.
  int bar_;
  // Work done so far.
  long work_done_;
  // Total work to do.
  long work_total_;
};

}  // namespace cs

#endif  // CS_PROGRESS_TABLE_H_
