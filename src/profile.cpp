/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "profile.h"

namespace cs
{

Profile::Profile()
{}

Profile::Profile(int ncols, int ndim) : RowMajorMatrix<float>(ncols, ndim)
{}

Profile::~Profile()
{}

void reset(Profile& profile)
{
    const int cols = profile.ncols();
    const int dim = profile.ndim();
    for(int i=0; i<cols; ++i)
        for(int j=0; i<dim; ++j)
            profile(i,j) = 0.0f;
}

}//cs
