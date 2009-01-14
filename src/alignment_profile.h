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
    // Constructs profile from serialized profile read from input stream.
    AlignmentProfile(std::istream& in, const SequenceAlphabet* alphabet);
    // Constructs a profile of the given sequence.
    explicit AlignmentProfile(const Sequence& sequence);
    // Constructs a profile of the given alignment with specified sequence weighting method.
    explicit AlignmentProfile(const Alignment& alignment, bool position_dependent_weights = true);
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

  protected:
    // Initializes the profile object with a serialized profile read from stream.
    virtual void unserialize(std::istream& in);
    // Prints the profile in serialization format to output stream.
    virtual void serialize(std::ostream& out) const;

  private:
    // Class identifier for serialization
    static const char kClass[];

    // Disallow copy and assign
    AlignmentProfile(const AlignmentProfile&);
    void operator=(const AlignmentProfile&);

    // Number of effective sequences in each alignment column.
    std::vector<float> neff_;
    // Flag indicating if the profile contains counts or (relative) frequencies.
    bool has_counts_;
};//AlignmentProfile

}//cs

#endif
