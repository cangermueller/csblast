// Copyright 2009, Andreas Biegert

#ifndef SRC_PROGRESS_TABLE_H_
#define SRC_PROGRESS_TABLE_H_

#include <cstdlib>
#include <cstdio>

namespace cs {

// Abstract base class for progress table to be subclassed for HMM training
// and EM clustering.
class ProgressTable {
 public:
  // To be used by derived classes.
  ProgressTable(FILE* fout = stdout, int width = 30);

  virtual ~ProgressTable() {}

  // Prints header information.
  virtual void print_header() = 0;
  // Starts a new row and prints statistics to outstream.
  virtual void print_row_begin() = 0;
  // Ends the current row and prints likelihood.
  virtual void print_row_end() = 0;
  // Sets total work per scan.
  void set_total_work(int work) { work_total_ = work; }
  // Resets the progress bar to zero.
  void reset();
  // Advances the progress bar proportional to the amount of work done.
  void print_progress(int work);

 protected:

  // Output stream.
  FILE* fout_;
  // Width of prograess bar.
  const int width_;
  // Currnt width of bar.
  int bar_;
  // Work done so far.
  int work_done_;
  // Total work to do.
  int work_total_;
};

}  // namespace cs

#endif  // SRC_PROGRESS_TABLE_H_
