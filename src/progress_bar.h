// Copyright 2009, Andreas Biegert

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
