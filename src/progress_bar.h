/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_PROGRESS_BAR_H_
#define CS_PROGRESS_BAR_H_

namespace cs {

// Simple graphiccal progress bar
class ProgressBar {
 public:
  // Constructs a progress bar of given length (incl. bracket delimiters)
  ProgressBar(FILE* fout, int width = 80, long todo = 0)
  : fout_(fout),
    steps_(width - 2),
    todo_(todo),
    prog_(0),
    done_(0) {}

  // Initializes total work todo to given value.
  void Init(long todo) { todo_ = todo; Reset(); }

  // Resets the progress bar to zero.
  void Reset() { done_ = 0; prog_ = 0; }

  // Advances the progress bar proportional to the amount of work done.
  void Advance(long work = 1) {
    assert_ne(0, todo_);
    int incr = static_cast<int>(round(static_cast<double>(done_ + work) / todo_ * steps_) - prog_);

    if (done_ == 0) fputc('[', fout_);
    fputs(std::string(incr, '=').c_str(), fout_);
    if (done_ + work == todo_) fputc(']', fout_);
    fflush(fout_);

    prog_ += incr;
    done_ += work;
  }

  void Complete() { while (done_ != todo_) Advance(); }

  // Returns output stream of the progress table
  FILE* stream() const { return fout_; }

 private:
  FILE* fout_;
  const int steps_;
  long todo_;
  int prog_;
  long done_;
};

}  // namespace cs

#endif  // CS_PROGRESS_BAR_H_
