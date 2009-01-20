#ifndef MATRIX_H
#define MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Matrix template class with continuous memory footprint in row major layout.

#include <vector>

#include "my_exception.h"

template<typename T>
class Matrix {
  public:
    Matrix();
    Matrix(int nrows, int ncols);
    Matrix(int nrows, int ncols, const T& val);

    // Access methods to get the (i,j) element:
    T&       operator() (int i, int j);
    const T& operator() (int i, int j) const;
    // Returns #rows in this matrix
    int nrows() const;
    // Returns #columns in this matrix
    int ncols() const;
    // Resizes the matrix to given dimensions. Old data is lost.
    void resize(int nrows, int ncols);

  private:
    int nrows_;
    int ncols_;
    std::vector<T> data_;
};



template<typename T>
Matrix<T>::Matrix()
        : nrows_(0), ncols_(0)
{}

template<typename T>
Matrix<T>::Matrix(int nrows, int ncols)
        : nrows_(nrows), ncols_(ncols), data_(nrows * ncols)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
}

template<typename T>
Matrix<T>::Matrix(int nrows, int ncols, const T& val)
        : nrows_(nrows), ncols_(ncols), data_(nrows * ncols, val)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
}

template<typename T>
inline int Matrix<T>::nrows() const
{ return nrows_; }

template<typename T>
inline int Matrix<T>::ncols() const
{ return ncols_; }

template<typename T>
inline T& Matrix<T>::operator() (int i, int j)
{ return data_[i*ncols_ + j]; }

template<typename T>
inline const T& Matrix<T>::operator() (int i, int j) const
{ return data_[i*ncols_ + j]; }

template<typename T>
void Matrix<T>::resize(int nrows, int ncols)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
    nrows_ = nrows;
    ncols_ = ncols;
    data_.resize(nrows * ncols);
}

#endif
