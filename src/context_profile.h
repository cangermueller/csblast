#ifndef CS_CONTEXT_PROFILE_H
#define CS_CONTEXT_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for profiles derived from alignments.

#include <iostream>
#include <vector>

#include "alignment.h"
#include "exception.h"
#include "profile.h"
#include "shared_ptr.h"
#include <vector>

namespace cs
{

template<class AlphabetType>
class ContextProfile : public Profile<AlphabetType>
{
  public:
    // Constructs profile from serialized profile read from input stream.
    ContextProfile(std::istream& in);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    ContextProfile(const Profile& other, int index, int length);

    virtual ~ContextProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<ContextProfile> > readall(std::istream& in);

    // Returns index of central profile column.
    int center() const { return (data_.size() - 1) / 2; }

  private:
    // Disallow copy and assign
    ContextProfile(const ContextProfile&);
    void operator=(const ContextProfile&);

    // Return serialization class identity.
    virtual const std::string class_identity() { static std::string id("ContextProfile"); return id;}
};  // ContextProfile



template<class AlphabetType>
ContextProfile<AlphabetType>::ContextProfile(std::istream& in, const SequenceAlphabet* alphabet)
        : Profile<AlphabetType>(alphabet)
{
    read(in);
    if (num_cols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but num_cols=%i!", num_cols());
}

template<class AlphabetType>
ContextProfile<AlphabetType>::ContextProfile(const Profile<AlphabetType>& other,
                                             int index,
                                             int length)
        : Profile<AlphabetType>(other, index, length)
{
    if (num_cols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but num_cols=%i!", num_cols());
}

template<class AlphabetType>
std::vector< shared_ptr< ContextProfile<AlphabetType> > > ContextProfile<AlphabetType>::readall(std::istream& in)
{
    std::vector< shared_ptr<ContextProfile> > profiles;
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<ContextProfile> p(new ContextProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

}  // cs

#endif
