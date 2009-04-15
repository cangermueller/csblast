// Copyright 2009, Andreas Biegert

#include "progress_table.h"

#include <cmath>
#include <cstdio>

#include <string>

namespace cs {

ProgressTable::ProgressTable(FILE* fout, int width)
    : fout_(fout),
      width_(width),
      bar_(0),
      work_done_(0),
      work_total_(0) {}

void ProgressTable::reset() {
  work_done_ = 0;
  bar_       = 0;
}

void ProgressTable::print_progress(int work) {
  if (work_total_ == 0) return;

  const int incr = round(static_cast<float>(work_done_ + work) /
                         work_total_ * (width_ - 2)) - bar_;
  if (work_done_ == 0) fputc('[', fout_);
  fputs(std::string(incr, '=').c_str(), fout_);
  if (work_done_ + work == work_total_) fputc(']', fout_);
  fflush(fout_);

  bar_       += incr;
  work_done_ += work;
}

};  // namespace cs
