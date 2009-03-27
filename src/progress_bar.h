#ifndef CS_PROGRESS_BAR_H
#define CS_PROGRESS_BAR_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Simple progress bar to be used in progress tables of HMM training and
// EM clustering.

#include <iostream>

namespace cs
{

class ProgressBar
{
  public:
    ProgressBar(std::ostream& out = std::cout, int width = 30)
            : width_(width),
              bar_(0),
              work_done_(0),
              work_total_(0),
              out_(out)
    { }

    // Sets the total amount of work to be done in each iteration.
    void set_total_work(int total)
    {
        work_total_ = total;
    }

    // Resets the progress bar to zero.
    void reset()
    {
        work_done_ = 0;
        bar_       = 0;
    }

    // Advances the progress bar proportional to the amount of work done.
    void print_progress(int work)
    {
        if (work_total_ == 0) return;

        const int incr = round(static_cast<float>(work_done_ + work) /
                               work_total_ * (width_ - 2)) - bar_;
        if (work_done_ == 0) out_ << "[";
        out_ << std::string(incr, '=');
        if (work_done_ + work == work_total_) out_ << "]";
        out_.flush();

        bar_       += incr;
        work_done_ += work;
    }

  protected:
    // Width of progress bar.
    const int width_;
    // Currnt width of bar.
    int bar_;
    // Work done so far.
    int work_done_;
    // Total work to do.
    int work_total_;
    // Output stream.
    std::ostream& out_;
};

}  // cs

#endif
