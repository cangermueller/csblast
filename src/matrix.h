// Copyright 2009, Andreas Biegert

#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

#include "stride_iter.h"

#include <valarray>
#include <numeric>
#include <algorithm>

using std::valarray;

template<class T>
class matrix {
 public:
  // public typedefs
  typedef T value_type;
  typedef matrix self;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  typedef value_type* row_type;
  typedef stride_iter<value_type*> col_type;
  typedef const value_type* const_row_type;
  typedef stride_iter<const value_type*> const_col_type;

  // constructors
  matrix() : num_rows_(0), num_cols_(0), m_() { }
  matrix(int r, int c) : num_rows_(r), num_cols_(c), m_(r * c) { }
  matrix(const self& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.m_) { }

  matrix(int r, int c, const T& val)
      : num_rows_(r), num_cols_(c), m_(r * c) {
    for (int i = 0; i < r * c; ++i) m_[i] = val;
  }

  template<typename New_T>
  explicit matrix(const valarray<New_T>& x)
      : num_rows_(x.size()), num_cols_(1), m_(x.size()) {
    for (int i =0 ; i < x.size(); ++i)
      m_[i] = x[i];
  }

  // allow construction from matricies of other types
  template<typename New_T>
  explicit matrix(const matrix<New_T>& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.size()) {
    copy(x.begin(), x.end(), m_.begin());
  }

  // public functions
  int num_rows() const { return num_rows_; }
  int num_cols() const { return num_cols_; }
  int size() const { return num_rows_ * num_cols_; }

  // element access
  row_type row_begin(int n) { return &m_[n * num_cols()]; }
  row_type row_end(int n) { return row_begin(n) + num_cols(); }
  col_type col_begin(int n) { return col_type(&m_[n], num_cols()); }
  col_type col_end(int n) { return col_begin(n) + num_cols(); }
  const_row_type row_begin(int n) const { return &m_[n * num_cols()]; }
  const_row_type row_end(int n) const { return row_begin(n) + num_cols(); }
  const_col_type col_begin(int n) const {
    return const_col_type(&m_[n], num_cols());
  }
  const_col_type col_end(int n) const { return col_begin(n) + num_cols(); }
  iterator begin() { return &m_[0]; }
  iterator end() { return begin() + size(); }
  const_iterator begin() const { return &m_[0]; }
  const_iterator end() const { return begin() + size(); }
  void resize(int r, int c, T v = T()) {
    m_.resize(r * c, v);
    num_rows_ = r;
    num_cols_ = c;
  }

  // operators
  self& operator=(const self& x) {
    m_.resize(x.size());
    m_ = x.m_;
    num_rows_ = x.num_rows_;
    num_cols_ = x.num_cols_;
    return *this;
  }
  self& operator=(value_type x) { m_ = x; return *this; }
  row_type operator[](int n) { return row_begin(n); }
  const_row_type operator[](int n) const { return row_begin(n); }
  self& operator+=(const self& x) { m_ += x.m_; return *this; }
  self& operator-=(const self& x) { m_ -= x.m_; return *this; }
  self& operator+=(value_type x) { m_ += x; return *this; }
  self& operator-=(value_type x) { m_ -= x; return *this; }
  self& operator*=(value_type x) { m_ *= x; return *this; }
  self& operator/=(value_type x) { m_ /= x; return *this; }
  self& operator%=(value_type x) { m_ %= x; return *this; }
  self operator-( ) { return -m_; }
  self operator+( ) { return +m_; }
  self operator!( ) { return !m_; }
  self operator~( ) { return ~m_; }

  // friend operators
  friend self operator+(const self& x, const self& y) { return self(x) += y; }
  friend self operator-(const self& x, const self& y) { return self(x) -= y; }
  friend self operator+(const self& x, value_type y) { return self(x) += y; }
  friend self operator-(const self& x, value_type y) { return self(x) -= y; }
  friend self operator*(const self& x, value_type y) { return self(x) *= y; }
  friend self operator/(const self& x, value_type y) { return self(x) /= y; }
  friend self operator%(const self& x, value_type y) { return self(x) %= y; }

 private:
  int num_rows_;
  int num_cols_;
  mutable valarray<T> m_;
};

#endif  // SRC_MATRIX_H_
