#ifndef CS_PROFILE_H
#define CS_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing plain columns of frequencies.

#include "row_major_matrix.h"

namespace cs
{

class Profile : private RowMajorMatrix<float>
{
public:
    Profile();
    Profile(int ncols, int ndim);
    virtual ~Profile();

    using RowMajorMatrix<float>::operator();
    using RowMajorMatrix<float>::resize;

    int ncols() const;
    int length() const;
    int ndim() const;
};//Profile



// Resets all entries in given profile to zero.
void reset(Profile& profile);

// Returns the number of columns in the profile.
inline int Profile::ncols() const
{ return RowMajorMatrix<float>::nrows(); }

// Returns the number of columns in the profile.
inline int Profile::length() const
{ return RowMajorMatrix<float>::nrows(); }

// Returns the number of entries per profile column.
inline int Profile::ndim() const
{ return RowMajorMatrix<float>::ncols(); }

}//cs

#endif
