/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "progress_table.h"

#include <cmath>

namespace cs
{

ProgressTable::ProgressTable(std::ostream& out, int width)
        : out_(out),
          width_(width),
          bar_(0),
          work_done_(0),
          work_total_(0)
{ }

void ProgressTable::reset()
{
    work_done_ = 0;
    bar_       = 0;
}

void ProgressTable::print_progress(int work)
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

};  // cs
