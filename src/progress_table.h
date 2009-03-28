#ifndef CS_PROGRESS_TABLE_H
#define CS_PROGRESS_TABLE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for progress table to be subclassed for HMM training
// and EM clustering.

#include <iostream>

namespace cs
{

class ProgressTable
{
  public:
     // To be used by derived classes.
    ProgressTable(std::ostream& out = std::cout,
                  int width = 30);

    virtual ~ProgressTable() { }

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
    // Returns reference to output stream.
    std::ostream& outstream() const { return out_; }

  protected:

    // Output stream.
    std::ostream& out_;
    // Width of prograess bar.
    const int width_;
    // Currnt width of bar.
    int bar_;
    // Work done so far.
    int work_done_;
    // Total work to do.
    int work_total_;
};

}  // cs

#endif
