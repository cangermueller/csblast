#ifndef STRIDE_ITER_H
#define STRIDE_ITER_H

#include <iterator>
#include <cassert>

template<class Iter_T>
class stride_iter
{
  public:
    // public typedefs
    typedef typename std::iterator_traits<Iter_T>::value_type value_type;
    typedef typename std::iterator_traits<Iter_T>::reference reference;
    typedef typename std::iterator_traits<Iter_T>::difference_type
    difference_type;
    typedef typename std::iterator_traits<Iter_T>::pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef stride_iter self;

    // constructors
    stride_iter( ) : m_(NULL), step_(0) { };
    stride_iter(const self& x) : m_(x.m_), step_(x.step_) { }
    stride_iter(Iter_T x, difference_type n) : m_(x), step_(n) { }

    // operators
    self& operator++( ) { m_ += step_; return *this; }
    self operator++(int) { self tmp = *this; m_ += step_; return tmp; }
    self& operator+=(difference_type x) { m_ += x * step_; return *this; }
    self& operator--( ) { m_ -= step_; return *this; }
    self operator--(int) { self tmp = *this; m_ -= step_; return tmp; }
    self& operator-=(difference_type x) { m_ -= x * step_; return *this; }
    reference operator[](difference_type n) { return m_[n * step_]; }
    reference operator*( ) { return *m_; }

    // friend operators
    friend bool operator==(const self& x, const self& y) {
        assert(x.step_ == y.step_);
        return x.m_ == y.m_;
    }
    friend bool operator!=(const self& x, const self& y) {
        assert(x.step_ == y.step_);
        return x.m_ != y.m_;
    }
    friend bool operator<(const self& x, const self& y) {
        assert(x.step_ == y.step_);
        return x.m_ < y.m_;
    }
    friend difference_type operator-(const self& x, const self& y) {
        assert(x.step_ == y.step_);
        return (x.m_ - y.m_) / x.step_;
    }
    friend self operator+(const self& x, difference_type y) {
        assert(x.step_ == y.step_);
        return x += y * x.step_;
    }
    friend self operator+(difference_type x, const self& y) {
        assert(x.step_ == y.step_);
        return y += x * x.step_;
    }
  private:
    Iter_T m_;
    difference_type step_;
};

#endif