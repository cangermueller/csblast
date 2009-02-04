#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "stride_iter.h"

#include <algorithm>
#include <google/sparsetable>
#include <numeric>

using google::sparsetable;

template<class T>
class sparse_matrix
{
  public:
    // public typedefs
    typedef T value_type;
    typedef sparse_matrix self;
    typedef typename sparsetable<value_type>::iterator iterator;
    typedef typename sparsetable<value_type>::const_iterator const_iterator;
    typedef typename sparsetable<value_type>::nonempty_iterator nonempty_iterator;
    typedef typename sparsetable<value_type>::const_nonempty_iterator const_nonempty_iterator;
    typedef typename sparsetable<value_type>::iterator row_type;
    typedef stride_iter< typename sparsetable<value_type>::iterator > col_type;
    typedef typename sparsetable<value_type>::const_iterator const_row_type;
    typedef stride_iter< typename sparsetable<value_type>::const_iterator > const_col_type;

    // constructors
    sparse_matrix() : num_rows_(0), num_cols_(0), m_() { }
    sparse_matrix(int r, int c) : num_rows_(r), num_cols_(c), m_(r * c) { }
    sparse_matrix(const self& x) : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.m_) { }

    // allow construction from matricies of other types
    template<typename New_T>
    explicit sparse_matrix(const sparse_matrix<New_T>& x)
            : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x)
    {}

    // public functions
    int num_rows() const { return num_rows_; }
    int num_cols() const { return num_cols_; }
    int size() const { return num_rows_ * num_cols_; }

    // element access
    row_type row_begin(int n) { return m_.begin() + n * num_cols(); }
    row_type row_end(int n) { return row_begin(n) + num_cols(); }
    col_type col_begin(int n) { return col_type(m_.begin() + n, num_cols()); }
    col_type col_end(int n) { return col_begin(n) + num_cols(); }
    const_row_type row_begin(int n) const { return m_.begin() + n * num_cols(); }
    const_row_type row_end(int n) const { return row_begin(n) + num_cols(); }
    const_col_type col_begin(int n) const { return col_type(m_.begin() + n, num_cols()); }
    const_col_type col_end(int n) const { return col_begin(n) + num_cols(); }
    iterator begin() { return m_.begin(); }
    iterator end() { return m_.end(); }
    const_iterator begin() const { return m_.begin(); }
    const_iterator end() const { return m_.end(); }
    nonempty_iterator nonempty_begin() { return m_.nonempty_begin(); }
    nonempty_iterator nonempty_end() { return m_.nonempty_end(); }
    const_nonempty_iterator nonempty_begin() const { return m_.nonempty_begin(); }
    const_nonempty_iterator nonempty_end() const { return m_.nonempty_end(); }
    void clear() { m_.clear(); }
    void resize(int r, int c)
    {
        m_.resize(r * c);
        num_rows_ = r;
        num_cols_ = c;
    }

    // operators
    self& operator=(const self& x)
    {
        m_.resize(x.size());
        m_ = x.m_;
        num_rows_ = x.num_rows_;
        num_cols_ = x.num_cols_;
        return *this;
    }
    row_type operator[](int n) { return row_begin(n); }
    const_row_type operator[](int n) const { return row_begin(n); }
    const T& get(int i, int j) { return m_.get(i*num_cols_ + j); }
    T& set(int i, int j, const T& val) { return m_.set(i*num_cols_ + j, val); }

  private:
    int num_rows_;
    int num_cols_;
    mutable sparsetable<T> m_;
};

#endif
