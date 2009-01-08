#ifndef CS_ROW_MAJOR_MATRIX_H
#define CS_ROW_MAJOR_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Matrix template class with continuous memory footprint in row major layout.

#include <vector>

#include "my_exception.h"

template<typename T>
class RowMajorMatrix {
public:
    RowMajorMatrix(int nrows, int ncols);

    // Access methods to get the (i,j) element:
    T&       operator() (int i, int j);
    const T& operator() (int i, int j) const;

    int nrows() const;  // #rows in this matrix
    int ncols() const;  // #columns in this matrix

private:
    int nrows_;
    int ncols_;
    std::vector<T> data_;
};

template<typename T>
inline int RowMajorMatrix<T>::nrows() const
{ return nrows_; }

template<typename T>
inline int RowMajorMatrix<T>::ncols() const
{ return ncols_; }

template<typename T>
inline T& RowMajorMatrix<T>::operator() (int row, int col)
{ return data_[row*ncols_ + col]; }

template<typename T>
inline const T& RowMajorMatrix<T>::operator() (int row, int col) const
{ return data_[row*ncols_ + col]; }

template<typename T>
RowMajorMatrix<T>::RowMajorMatrix(int nrows, int ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows * ncols)
{
    if (nrows == 0 || ncols == 0)
        throw MyException("Bad size arguments for matrix: nrows=%i ncols=%i", nrows, ncols);
}

#endif
