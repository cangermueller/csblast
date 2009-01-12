#ifndef CS_SEQUENCE_PROFILE_H
#define CS_SEQUENCE_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies over a sequence alphabet.

#include <iostream>
#include <vector>

#include "smart_ptr.h"

namespace cs
{

// Forward declarations
class Sequence;
class SequenceAlphabet;

class SequenceProfile
{
  public:
    friend std::istream& operator>> (std::istream& in, SequenceProfile& profile);
    friend std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile);

    // Constructs a profile with ncols columns over alphabet with entries initialized to zero.
    SequenceProfile(int ncols, const SequenceAlphabet* alphabet);
    // Constructs profile from serialized profile read from input stream.
    SequenceProfile(std::istream& in, const SequenceAlphabet* alphabet);
    // Constructs a profile of the given sequence.
    explicit SequenceProfile(const Sequence& sequence);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    SequenceProfile(const SequenceProfile& other, int index, int length);

    virtual ~SequenceProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< SmartPtr<SequenceProfile> > read(std::istream& in,
                                                         const SequenceAlphabet* alphabet);
    // Access methods to get the (i,j) element
    float&       operator() (int i, int j) { return data_[i*ndim_ + j]; }
    const float& operator() (int i, int j) const { return data_[i*ndim_ + j]; }
    // Returns #rows in this matrix
    int ncols() const { return ncols_; }
    // Returns #columns in this matrix
    int ndim() const { return ndim_; }
    // Returns the underlying sequence alphabet of the profile.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

  private:
    // Scaling factor for serialization of profile log values
    static const int kScaleFactor = 1000;

    // Disallow copy and assign
    SequenceProfile(const SequenceProfile&);
    void operator=(const SequenceProfile&);

    // Initializes the profile object with a serialized profile read from stream.
    void init(std::istream& in);
    // Resize the profile matrix to given dimensions. Attention: old data is lost!
    void resize(int ncols, int ndim);

    // Number of columns in the profile
    int ncols_;
    // Number of entries per column
    int ndim_;
    // Profile matrix in row major format
    std::vector<float> data_;
    // Sequence alphabet over which the profile records probabilities. Note that the profile
    // does not include the 'any' character (ndim = alphabet_.size()-1).
    const SequenceAlphabet* alphabet_;
};//SequenceProfile



// Initializes a sequence profile from seralized profile in input stream.
std::istream& operator>> (std::istream& in, SequenceProfile& profile);

// Prints a sequence profile in human readable serialization format.
std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile);

// Resets all entries in given profile to the given value or zero if none is given.
void reset(SequenceProfile& profile, float value = 0.0f);

}//cs

#endif
