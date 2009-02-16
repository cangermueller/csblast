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

namespace cs
{

template<class Alphabet_T>
class ContextProfile : public Profile<Alphabet_T>
{
  public:
    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::num_cols;

     // Constructs a dummy context profile.
    ContextProfile() {}
    // Constructs a context profile with num_cols columns initialized to zero.
    explicit ContextProfile(int num_cols);
    // Constructs profile from serialized profile read from input stream.
    explicit ContextProfile(std::istream& in);
    // Constructs a context profile from simple profile and checks if length is valid.
    explicit ContextProfile(const Profile<Alphabet_T>& profile);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    ContextProfile(const Profile<Alphabet_T>& other, int index, int length);

    virtual ~ContextProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<ContextProfile> > readall(std::istream& in);

    // Returns index of central profile column.
    int center() const { return (num_cols() - 1) / 2; }

  private:
    // Disallow copy and assign
    ContextProfile(const ContextProfile&);
    void operator=(const ContextProfile&);

    // Checks if profile has odd number of columns.
    void check();

    // Return serialization class identity.
    virtual const std::string class_identity() { static std::string id("ContextProfile"); return id;}
};  // ContextProfile



template<class Alphabet_T>
ContextProfile<Alphabet_T>::ContextProfile(int num_cols)
        : Profile<Alphabet_T>(num_cols)
{
    check();
}

template<class Alphabet_T>
ContextProfile<Alphabet_T>::ContextProfile(const Profile<Alphabet_T>& profile)
        : Profile<Alphabet_T>(profile)
{
    check();
}

template<class Alphabet_T>
ContextProfile<Alphabet_T>::ContextProfile(std::istream& in)
        : Profile<Alphabet_T>(in)
{
    check();
}

template<class Alphabet_T>
ContextProfile<Alphabet_T>::ContextProfile(const Profile<Alphabet_T>& other,
                                           int index,
                                           int length)
        : Profile<Alphabet_T>(other, index, length)
{
    check();
}

template<class Alphabet_T>
std::vector< shared_ptr< ContextProfile<Alphabet_T> > > ContextProfile<Alphabet_T>::readall(std::istream& in)
{
    std::vector< shared_ptr<ContextProfile> > profiles;
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<ContextProfile> p(new ContextProfile(in));
        profiles.push_back(p);
    }

    return profiles;
}

template<class Alphabet_T>
void ContextProfile<Alphabet_T>::check()
{
    if (num_cols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but num_cols=%i!", num_cols());
}

}  // cs

#endif
