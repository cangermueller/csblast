// Copyright 2009, Andreas Biegert

#ifndef SRC_COUNT_PROFILE_H_
#define SRC_COUNT_PROFILE_H_

#include <cmath>
#include <cstdio>
#include <cstring>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "alignment.h"
#include "exception.h"
#include "matrix.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs {

// A container class for profiles derived from alignments.
template<class Alphabet>
class CountProfile : public Profile<Alphabet> {
 public:
  // Needed to access names in templatized Profile base class
  using Profile<Alphabet>::num_cols;
  using Profile<Alphabet>::alphabet_size;
  using Profile<Alphabet>::read;
  using Profile<Alphabet>::logspace;
  using Profile<Alphabet>::transform_to_logspace;
  using Profile<Alphabet>::transform_to_linspace;

  // Constructs profile from serialized profile read from input stream.
  explicit CountProfile(std::istream& in);
  // Constructs profile from serialized profile read from input stream.
  explicit CountProfile(FILE* fin);
  // Constructs a profile of the given sequence.
  explicit CountProfile(const Sequence<Alphabet>& sequence);
  // Constructs a profile of the given alignment with specified sequence weighting method.
  explicit CountProfile(const Alignment<Alphabet>& alignment,
                        bool position_specific_weights = true);
  // Creates a profile from subprofile starting at column index and length columns long.
  CountProfile(const CountProfile& other, int index, int length);

  virtual ~CountProfile() { }

  // Reads all available profiles from the input stream and returns them in a vector.
  static void readall(std::istream& in, std::vector< shared_ptr<CountProfile> >& v);
  // Reads all available profiles from the input stream and returns them in a vector.
  static void readall(FILE* in, std::vector< shared_ptr<CountProfile> >* v);
  // Returns the number of effective sequences in alignment column i
  float neff(int i) const { return neff_[i]; }
  // Converts the profile to counts of alphabet letters.
  void convert_to_counts();
  // Converts the profile back to relative frequencies of alphabet letters.
  void convert_to_frequencies();
  // Returns true if the profile contains counts.
  bool has_counts() const { return has_counts_; }

 protected:
  // Needed to access names in templatized Profile base class
  using Profile<Alphabet>::data_;
  using Profile<Alphabet>::SCALE_FACTOR;

  using Profile<Alphabet>::LOG_SCALE;
  using Profile<Alphabet>::BUFFER_SIZE;

  // Reads and initializes serialized scalar data members from stream.
  virtual void read_header(std::istream& in);
  // Reads and initializes serialized scalar data members from stream.
  virtual void read_header(FILE* fin);
  // Reads and initializes array data members from stream.
  virtual void read_body(std::istream& in);
  // Reads and initializes array data members from stream.
  virtual void read_body(FILE* fin);
  // Writes serialized scalar data members to stream.
  virtual void write_header(std::ostream& out) const;
  // Writes serialized array data members to stream.
  virtual void write_body(std::ostream& out) const;
  // Prints the profile in human-readable format to output stream.
  virtual void print(std::ostream& out) const;

 private:
  // Class identifier
  static const char* CLASS_ID;

  // Return serialization class identity.
  virtual const std::string class_identity() const {
    static std::string id("CountProfile");
    return id;
  }
  virtual const char* class_id() const { return CLASS_ID; }

  // Number of effective sequences in each alignment column.
  std::vector<float> neff_;
  // Flag indicating if the profile contains counts or (relative) frequencies.
  bool has_counts_;
};  // CountProfile

}  // cs

#endif  // SRC_COUNT_PROFILE_H_
