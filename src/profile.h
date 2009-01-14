#ifndef CS_PROFILE_H
#define CS_PROFILE_H
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

class Profile
{
  public:
    friend std::istream& operator>> (std::istream& in, Profile& profile);
    friend std::ostream& operator<< (std::ostream& out, const Profile& profile);

    // Constructs a dummy profile with given alphabet.
    Profile(const SequenceAlphabet* alphabet);
    // Constructs a profile with ncols columns over alphabet with entries initialized to zero.
    Profile(int ncols, const SequenceAlphabet* alphabet);
    // Constructs profile from serialized profile read from input stream.
    Profile(std::istream& in, const SequenceAlphabet* alphabet);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    Profile(const Profile& other, int index, int length);

    virtual ~Profile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< SmartPtr<Profile> > read(std::istream& in,
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

  protected:
    // Scaling factor for serialization of profile log values
    static const int kScaleFactor = 1000;

    // Initializes the profile object with a serialized profile read from stream.
    virtual void unserialize(std::istream& in);
    // Prints the profile in serialization format to output stream.
    virtual void serialize(std::ostream& out) const;
    // Resize the profile matrix to given dimensions. Attention: old data is lost!
    void resize(int ncols, int ndim);

  private:
    // Class identifier for serialization
    static const char kClass[];

    // Disallow copy and assign
    Profile(const Profile&);
    void operator=(const Profile&);

    // Number of columns in the profile
    int ncols_;
    // Number of entries per column
    int ndim_;
    // Profile matrix in row major format
    std::vector<float> data_;
    // Sequence alphabet over which the profile records probabilities. Note that the profile
    // does not include the 'any' character (ndim = alphabet_.size()-1).
    const SequenceAlphabet* alphabet_;
};//Profile



// Initializes a sequence profile from seralized profile in input stream.
std::istream& operator>> (std::istream& in, Profile& profile);

// Prints a sequence profile in human readable serialization format.
std::ostream& operator<< (std::ostream& out, const Profile& profile);

// Resets all entries in given profile to the given value or zero if none is given.
void reset(Profile& profile, float value = 0.0f);

// Normalize profile columns to value or to one if none provided.
void normalize(Profile& profile, float value = 1.0f);

}//cs

#endif
