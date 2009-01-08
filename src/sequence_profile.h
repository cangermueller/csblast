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
    friend std::istream& operator>> (std::istream& i, SequenceProfile& profile);
    friend std::ostream& operator<< (std::ostream& o, const SequenceProfile& profile);

    SequenceProfile(int ncols,
                    const SequenceAlphabet* alphabet);
    virtual ~SequenceProfile();

    using RowMajorMatrix<float>::operator();
    int ncols() const;
    int length() const;
    int nalph() const;
    const SequenceAlphabet& alphabet() const;

private:
    const SequenceAlphabet* alphabet_;
};//SequenceProfile



std::istream& operator>> (std::istream& i, SequenceProfile& profile);

std::ostream& operator<< (std::ostream& o, const SequenceProfile& profile);

void normalize(SequenceProfile& profile);


inline int SequenceProfile::ncols() const
{ return RowMajorMatrix<float>::nrows(); }

inline int SequenceProfile::length() const
{ return RowMajorMatrix<float>::nrows(); }

inline int SequenceProfile::nalph() const
{ return RowMajorMatrix<float>::ncols(); }

inline const SequenceAlphabet& SequenceProfile::alphabet() const
{ return *alphabet_; }

}//cs

#endif
