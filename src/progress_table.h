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

// Forward declaration;
template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ExpectationMaximization;


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ProgressTable
{
  public:
    typedef ExpectationMaximization<Alphabet_T, Subject_T> em_type;

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
    // To be used by derived classes.
    ProgressTable(const em_type* em,
                  std::ostream& out = std::cout,
                  int width = 30);

    // Pointer to EM object.
    const em_type* em_;
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



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
ProgressTable<Alphabet_T, Subject_T>::ProgressTable(const em_type* em,
                                                    std::ostream& out,
                                                    int width)
        : em_(em),
          out_(out),
          width_(width),
          bar_(0),
          work_done_(0),
          work_total_(0)
{ }

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ProgressTable<Alphabet_T, Subject_T>::reset()
{
    work_done_ = 0;
    bar_       = 0;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ProgressTable<Alphabet_T, Subject_T>::print_progress(int work)
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

}  // cs

#endif
