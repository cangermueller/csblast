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

#include "matrix.h"
#include "shared_ptr.h"

namespace cs
{

// Forward declarations
class Sequence;
class SequenceAlphabet;

class Profile
{
  public:
    typedef matrix<float>::row_type col_type;
    typedef matrix<float>::const_row_type const_col_type;
    typedef matrix<float>::iterator iterator;
    typedef matrix<float>::const_iterator const_iterator;

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
    static std::vector< shared_ptr<Profile> > readall(std::istream& in,
                                                      const SequenceAlphabet* alphabet);

    // Access methods to get the (i,j) element
    col_type operator[](int i) { return data_[i]; }
    const_col_type operator[](int i) const { return data_[i]; }
    // Returns #columns in the profile
    int ncols() const { return data_.nrows(); }
    // Returns #entries per column
    int nalph() const { return data_.ncols(); }
    // Returns the total number of elements in the profile.
    int size() const { return data_.size(); }
    // Transforms profile to logspace
    virtual void transform_to_logspace();
    // Transforms profile to linspace
    virtual void transform_to_linspace();
    // Returns true if the profile is in logspace
    bool logspace() const { return logspace_; }
    // Returns the underlying sequence alphabet of the profile.
    const SequenceAlphabet* alphabet() const { return alphabet_; }
    // Returns an iterator to the first element in profile column i.
    col_type col_begin(int i) { return data_.row_begin(i); }
    // Returns an iterator just past the end of profile column i.
    col_type col_end(int i) { return data_.row_end(i); }
    // Returns a const iterator to the first element in profile column i.
    const_col_type col_begin(int i) const { return data_.row_begin(i); }
    // Returns a const iterator just past the end of profile column i.
    const_col_type col_end(int i) const { return data_.row_end(i); }
    // Returns an iterator to the first element in the profile matrix.
    iterator begin() { return data_.begin(); }
    // Returns an iterator just past the end of the profile matrix.
    iterator end() { return data_.end(); }
    // Returns a const iterator to the first element in the profile matrix.
    const_iterator begin() const { return data_.begin(); }
    // Returns a const iterator just past the end of the profile matrix.
    const_iterator end() const { return data_.end(); }
    // Initializes the profile object with a serialized profile read from stream.
    void read(std::istream& in);
    // Writes the profile in serialization format to output stream.
    void write(std::ostream& out) const;

    friend std::ostream& operator<< (std::ostream& out, const Profile& p);

  protected:
    // Scaling factor for serialization of profile log values
    static const int kScaleFactor = 1000;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized array data members to stream.
    virtual void write_body(std::ostream& out) const;
    // Resize the profile matrix to given dimensions. Attention: old data is lost!
    void resize(int ncols, int nalph);
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

     // Profile matrix in row major format
    matrix<float> data_;
    // Flag indicating if profile is in log- or linspace
    bool logspace_;
    // Sequence alphabet over which the profile records probabilities are recorded.
    const SequenceAlphabet* alphabet_;

  private:
    // Disallow copy and assign
    Profile(const Profile&);
    void operator=(const Profile&);

    // Returns serialization class identity.
    virtual const std::string& class_identity() const { static std::string id("Profile"); return id;}
};  // Profile



// Resets all entries in given profile to the given value or zero if none is given.
void reset(Profile& profile, float value = 0.0f);

// Normalize profile columns to value or to one if none provided.
void normalize(Profile& profile, float value = 1.0f);

// Prints profile in human-readable format for debugging.
inline std::ostream& operator<< (std::ostream& out, const Profile& p)
{
    p.print(out);
    return out;
}

}  // cs

#endif
