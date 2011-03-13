// Copyright 2009, Andreas Biegert

#ifndef CS_SPARSE_MATRIX_H_
#define CS_SPARSE_MATRIX_H_

#include "stride_iter.h"

#include <algorithm>
#include <google/sparsetable>
#include <numeric>

using google::sparsetable;

template<class T>
class SparseMatrix {
 public:
  // Public typedefs
  typedef T value_type;
  typedef SparseMatrix self;
  typedef typename sparsetable<T>::iterator iterator;
  typedef typename sparsetable<T>::const_iterator const_iterator;
  typedef typename sparsetable<T>::nonempty_iterator nonempty_iterator;
  typedef typename sparsetable<T>::const_nonempty_iterator const_nonempty_iterator;
  typedef typename sparsetable<T>::iterator row_type;
  typedef stride_iter< typename sparsetable<T>::iterator > col_type;
  typedef typename sparsetable<T>::const_iterator const_row_type;
  typedef stride_iter< typename sparsetable<T>::const_iterator > const_col_type;

  // Constructors
  SparseMatrix() : num_rows_(0), num_cols_(0), m_() {}
  SparseMatrix(int r, int c) : num_rows_(r), num_cols_(c), m_(r * c) {}
  SparseMatrix(const self& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.m_) {}

  // Allow construction from matricies of other types
  template<typename New_T>
  explicit SparseMatrix(const SparseMatrix<New_T>& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x) {}

  // public functions
  int num_rows() const { return num_rows_; }
  int num_cols() const { return num_cols_; }
  int size() const { return num_rows_ * num_cols_; }
  int num_nonempty() const { return m_.num_nonempty(); }

  // Element access
  row_type row_begin(int n) { return m_.begin() + n * num_cols(); }
  row_type row_end(int n) { return row_begin(n) + num_cols(); }
  col_type col_begin(int n) { return col_type(m_.begin() + n, num_cols()); }
  col_type col_end(int n) { return col_begin(n) + num_cols(); }
  const_row_type row_begin(int n) const { return m_.begin() + n * num_cols(); }
  const_row_type row_end(int n) const { return row_begin(n) + num_cols(); }
  const_col_type col_begin(int n) const {
    return col_type(m_.begin() + n, num_cols());
  }
  const_col_type col_end(int n) const { return col_begin(n) + num_cols(); }
  iterator begin() { return m_.begin(); }
  iterator end() { return m_.end(); }
  const_iterator begin() const { return m_.begin(); }
  const_iterator end() const { return m_.end(); }
  nonempty_iterator nonempty_begin() { return m_.nonempty_begin(); }
  nonempty_iterator nonempty_end() { return m_.nonempty_end(); }
  const_nonempty_iterator nonempty_begin() const { return m_.nonempty_begin(); }
  const_nonempty_iterator nonempty_end() const { return m_.nonempty_end(); }
  void erase(int r, int c) { m_.erase(r*num_cols_ + c); }
  void clear() { m_.clear(); }
  const T& get(int r, int c) const { return m_.get(r*num_cols_ + c); }
  T& set(int r, int c, const T& val) { return m_.set(r*num_cols_ + c, val); }
  bool test(int r, int c) const { return m_.test(r*num_cols_ + c); }
  void Resize(int r, int c) {
    m_.resize(r * c);
    num_rows_ = r;
    num_cols_ = c;
  }

  // Operators
  self& operator=(const self& x) {
    m_.resize(x.size());
    m_ = x.m_;
    num_rows_ = x.num_rows_;
    num_cols_ = x.num_cols_;
    return *this;
  }
  row_type operator[](int n) { return row_begin(n); }
  const_row_type operator[](int n) const { return row_begin(n); }

 private:
  int num_rows_;
  int num_cols_;
  mutable sparsetable<T> m_;
};

#endif  // CS_SPARSE_MATRIX_H_
