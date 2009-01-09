#ifndef CS_SEQUENCE_PROFILE_H
#define CS_SEQUENCE_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies over a sequence alphabet.

#include <iostream>
#include <streambuf>

#include "profile.h"
#include "sequence.h"
#include "sequence_alphabet.h"

namespace cs
{

class SequenceProfile : public Profile
{
public:
    friend std::istream& operator>> (std::istream& in, SequenceProfile& profile);

    SequenceProfile(int ncols, int ndim);
    SequenceProfile(int ncols,
                    const SequenceAlphabet* alphabet);
    SequenceProfile(const Sequence& sequence);
    virtual ~SequenceProfile();

    const SequenceAlphabet& alphabet() const;

private:
    const SequenceAlphabet* alphabet_;
};//SequenceProfile


// Initializes a sequence profile from seralized profile in input stream.
std::istream& operator>> (std::istream& in, SequenceProfile& profile);

// Prints a sequence profile in human readable serialization format.
std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile);


// Returns a reference to the underlying alphabet of the profile.
inline const SequenceAlphabet& SequenceProfile::alphabet() const
{ return *alphabet_; }

}//cs

#endif
