#ifndef CS_ALIGNMENT_PROFILE_H
#define CS_ALIGNMENT_PROFILE_H
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
#include "smart_ptr.h"

namespace cs
{

// Forward declarations
class SequenceAlphabet;

class AlignmentProfile : public Profile
{
  public:
    friend std::istream& operator>> (std::istream& in, AlignmentProfile& profile);
    friend std::ostream& operator<< (std::ostream& out, const AlignmentProfile& profile);

    // Constructs profile from serialized profile read from input stream.
    AlignmentProfile(std::istream& in, const SequenceAlphabet* alphabet);
    // Constructs a profile of the given alignment.
    explicit AlignmentProfile(const Alignment& alignment);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    AlignmentProfile(const AlignmentProfile& other, int index, int length);

    virtual ~AlignmentProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< SmartPtr<AlignmentProfile> > read(std::istream& in,
                                                          const SequenceAlphabet* alphabet);
    // Returns the number of effective sequences in alignment column i
    float neff(int i) const { return neff_[i]; }
    // Converts the profile to counts of alphabet letters.
    void convert_to_counts();
    // Converts the profile back to relative frequencies of alphabet letters.
    void convert_to_frequencies();
    // Returns true if the profile contains counts.
    bool has_counts() const { return has_counts_; }

  private:
    // Disallow copy and assign
    AlignmentProfile(const AlignmentProfile&);
    void operator=(const AlignmentProfile&);

    // Initializes the profile object with a serialized profile read from stream.
    virtual void init(std::istream& in);

    // Number of effective sequences in each alignment column.
    std::vector<float> neff_;
    // Flag indicating if the profile contains counts or (relative) frequencies.
    bool has_counts_;
};//AlignmentProfile



// Initializes an alignment profile from seralized profile in input stream.
std::istream& operator>> (std::istream& in, AlignmentProfile& profile);

// Prints an alignment profile in human readable serialization format.
std::ostream& operator<< (std::ostream& out, const AlignmentProfile& profile);

}//cs

#endif
