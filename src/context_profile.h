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
#include "profile.h"
#include "shared_ptr.h"

namespace cs
{

// Forward declarations
class SequenceAlphabet;

class ContextProfile : public Profile
{
  public:
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    ContextProfile(const Profile& other, int index, int length);
    // Constructs profile from serialized profile read from input stream.
    ContextProfile(std::istream& in, const SequenceAlphabet* alphabet);

    virtual ~ContextProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<ContextProfile> > readall(std::istream& in,
                                                             const SequenceAlphabet* alphabet);

    // Returns index of central profile column.
    int center() const { return (data_.size() - 1) / 2; }

  private:
    // Disallow copy and assign
    ContextProfile(const ContextProfile&);
    void operator=(const ContextProfile&);

    // Return serialization class identity.
    virtual const std::string& class_identity() { static std::string id("ContextProfile"); return id;}
};  // ContextProfile

}  // cs

#endif
