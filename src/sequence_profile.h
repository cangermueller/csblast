#ifndef CS_SEQUENCE_PROFILE_H
#define CS_SEQUENCE_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies for sequence alphabet

#include "row_major_matrix.h"
#include "sequence_alphabet.h"

namespace cs
{

class SequenceProfile : private RowMajorMatrix<float>
{
public:
    SequenceProfile(size_t ncols, const SequenceAlphabet* const alphabet);
    virtual ~SequenceProfile();

    using RowMajorMatrix<float>::operator();
    size_t ncols() const;
    size_t length() const;
    size_t nalph() const;

    const SequenceAlphabet& alphabet();

private:
    const SequenceAlphabet* const alphabet_;
};

inline size_t SequenceProfile::ncols() const
{ return RowMajorMatrix<float>::nrows(); }

inline size_t SequenceProfile::length() const
{ return RowMajorMatrix<float>::nrows(); }

inline size_t SequenceProfile::nalph() const
{ return RowMajorMatrix<float>::ncols(); }

void normalize(SequenceProfile& profile);

}//cs

#endif
