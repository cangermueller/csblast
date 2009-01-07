#ifndef CS_SEQUENCE_PROFILE_H
#define CS_SEQUENCE_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies for sequence alphabet

#include "Row_major_matrix.h"
#include "Sequence_alphabet.h"

class Sequence_profile : private Row_major_matrix<float>
{
public:
    Sequence_profile(size_t ncols, const Sequence_alphabet* const alphabet);
    virtual ~Sequence_profile();

    using Row_major_matrix<float>::operator();
    size_t ncols() const;
    size_t length() const;
    size_t nalph() const;

    const Sequence_alphabet& alphabet();

private:
    const Sequence_alphabet* const alphabet_;
};

inline size_t Sequence_profile::ncols() const
{ return Row_major_matrix<float>::nrows(); }

inline size_t Sequence_profile::length() const
{ return Row_major_matrix<float>::nrows(); }

inline size_t Sequence_profile::nalph() const
{ return Row_major_matrix<float>::ncols(); }

void normalize(Sequence_profile& profile);

#endif
