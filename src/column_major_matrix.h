#ifndef COLUMN_MAJOR_MATRIX_H
#define COLUMN_MAJOR_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Matrix template class with continuous memory footprint in column major layout.

#include <vector>

#include "my_exception.h"

template<typename T>
class ColumnMajorMatrix {
  public:
    ColumnMajorMatrix(int nrows, int ncols);

    // Access methods to get the (i,j) element:
    T&       operator() (int i, int j);
    const T& operator() (int i, int j) const;

    int nrows() const;  // #rows in this matrix
    int ncols() const;  // #columns in this matrix

    // Resize the matrix to given dimensions. Old data is NOT retained.
    void resize(int nrows, int ncols);

  protected:
    ColumnMajorMatrix();

  private:
    int nrows_;
    int ncols_;
    std::vector<T> data_;
};


template<typename T>
ColumnMajorMatrix<T>::ColumnMajorMatrix()
    : nrows_(0), ncols_(0), data_(0)
{}

template<typename T>
ColumnMajorMatrix<T>::ColumnMajorMatrix(int nrows, int ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows * ncols)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
}

template<typename T>
inline int ColumnMajorMatrix<T>::nrows() const
{ return nrows_; }

template<typename T>
inline int ColumnMajorMatrix<T>::ncols() const
{ return ncols_; }

template<typename T>
inline T& ColumnMajorMatrix<T>::operator() (int row, int col)
{ return data_[row + col*nrows_]; }

template<typename T>
inline const T& ColumnMajorMatrix<T>::operator() (int row, int col) const
{ return data_[row + col*nrows_]; }

template<typename T>
void ColumnMajorMatrix<T>::resize(int nrows, int ncols)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
    nrows_ = nrows;
    ncols_ = ncols;
    data_.resize(nrows * ncols);
}

#endif
