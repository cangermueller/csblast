/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "context_profile.h"

#include <vector>

#include "exception.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs
{

ContextProfile::ContextProfile(std::istream& in, const SequenceAlphabet* alphabet)
        : Profile(alphabet)
{
    read(in);
    if (ncols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but ncols=%i!", ncols());
}

ContextProfile::ContextProfile(const Profile& other,
                              int index,
                              int length)
        : Profile(other, index, length)
{
    if (ncols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but ncols=%i!", ncols());
}

std::vector< shared_ptr<ContextProfile> > ContextProfile::readall(std::istream& in,
                                                                 const SequenceAlphabet* alphabet)
{
    std::vector< shared_ptr<ContextProfile> > profiles;
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<ContextProfile> p(new ContextProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

}  // cs
