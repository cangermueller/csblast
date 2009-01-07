#ifndef CS_ROW_MAJOR_MATRIX_H
#define CS_ROW_MAJOR_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Matrix template class with continuous memory footprint in row major layout.

#include <cstddef>
#include <vector>

template<typename T>
class Matrix {
public:
    Matrix(size_t nrows, size_t ncols);

    // Access methods to get the (i,j) element:
    T&       operator() (size_t i, size_t j);
    const T& operator() (size_t i, size_t j) const;

    size_t nrows() const;  // #rows in this matrix
    size_t ncols() const;  // #columns in this matrix

private:
    size_t nrows_;
    size_t ncols_;
    std::vector<T> data_;
};

template<typename T>
inline size_t Matrix<T>::nrows() const
{ return nrows_; }

template<typename T>
inline size_t Matrix<T>::ncols() const
{ return ncols_; }

template<typename T>
inline T& Matrix<T>::operator() (size_t row, size_t col)
{ return data_[row*ncols_ + col]; }

template<typename T>
inline const T& Matrix<T>::operator() (size_t row, size_t col) const
{ return data_[row*ncols_ + col]; }

template<typename T>
Matrix<T>::Matrix(size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows * ncols)
{ }

#endif
