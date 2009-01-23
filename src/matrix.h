#ifndef MATRIX_H
#define MATRIX_H

#include "stride_iter.h"

#include <valarray>
#include <numeric>
#include <algorithm>

using std::valarray;

template<class Value_T>
class matrix
{
  public:
    // public typedefs
    typedef Value_T value_type;
    typedef matrix self;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef Value_T* row_type;
    typedef stride_iter<value_type*> col_type;
    typedef const value_type* const_row_type;
    typedef stride_iter<const value_type*> const_col_type;

    // constructors
    matrix() : nrows_(0), ncols_(0), m_() { }
    matrix(int r, int c) : nrows_(r), ncols_(c), m_(r * c) { }
    matrix(const self& x) : nrows_(x.nrows_), ncols_(x.ncols_), m_(x.m_) { }

    matrix(int r, int c, const Value_T& val)
            : nrows_(r), ncols_(c), m_(r * c)
    {
        for (int i = 0; i < r * c; ++i) m_[i] = val;
    }

    template<typename T>
    explicit matrix(const valarray<T>& x)
            : nrows_(x.size()), ncols_(1), m_(x.size() + 1)
    {
        for (int i =0 ; i < x.size(); ++i) m_[i] = x[i];
    }

    // allow construction from matricies of other types
    template<typename T>
    explicit matrix(const matrix<T>& x)
            : nrows_(x.nrows_), ncols_(x.ncols_), m_(x.size() + 1)
    {
        copy(x.begin(), x.end(), m_.begin());
    }

    // public functions
    int nrows() const { return nrows_; }
    int ncols() const { return ncols_; }
    int size() const { return nrows_ * ncols_; }

    // element access
    row_type row_begin(int n) { return &m_[n * ncols()]; }
    row_type row_end(int n) { return row_begin(n) + ncols(); }
    col_type col_begin(int n) { return col_type(&m_[n], ncols()); }
    col_type col_end(int n) { return col_begin(n) + ncols(); }
    const_row_type row_begin(int n) const { return &m_[n * ncols()]; }
    const_row_type row_end(int n) const { return row_begin(n) + ncols( ); }
    const_col_type col_begin(int n) const { return col_type(&m_[n], ncols( )); }
    const_col_type col_end(int n) const { return col_begin(n) + ncols( ); }
    iterator begin() { return &m_[0]; }
    iterator end() { return begin() + size(); }
    const_iterator begin() const { return &m_[0]; }
    const_iterator end() const { return begin() + size( ); }
    void resize(int r, int c)
    {
        m_.resize(r * c);
        nrows_ = r;
        ncols_ = c;
    }

    // operators
    self& operator=(const self& x)
    {
        m_.resize(x.size());
        m_ = x.m_;
        nrows_ = x.nrows_;
        ncols_ = x.ncols_;
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
    int nrows_;
    int ncols_;
    mutable valarray<Value_T> m_;
};

#endif
