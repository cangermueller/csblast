// Copyright 2009, Andreas Biegert

#ifndef SRC_COUNT_PROFILE_H_
#define SRC_COUNT_PROFILE_H_

#include <cmath>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "alignment.h"
#include "matrix.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"

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
  explicit CountProfile(FILE* fin);
  // Constructs a profile of the given sequence.
  explicit CountProfile(const Sequence<Alphabet>& sequence);
  // Constructs a profile of the given alignment with specified sequence
  // weighting method.
  explicit CountProfile(const Alignment<Alphabet>& alignment,
                        bool position_specific_weights = true);
  // Creates a profile from subprofile starting at column index and length
  // columns long.
  CountProfile(const CountProfile& other, int index, int length);

  virtual ~CountProfile() {}

  // Reads all available profiles from the input stream and returns them in a
  // vector.
  static void readall(FILE* in, std::vector< shared_ptr<CountProfile> >* v);
  // Returns the number of effective sequences in alignment column i
  float neff(int i) const { return neff_[i]; }

 protected:
  // Needed to access names in templatized Profile base class
  using Profile<Alphabet>::data_;
  using Profile<Alphabet>::kLogScale;
  using Profile<Alphabet>::kBufferSize;

  // Reads and initializes serialized scalar data members from stream.
  virtual void read_header(FILE* fin);
  // Reads and initializes array data members from stream.
  virtual void read_body(FILE* fin);
  // Writes serialized array data members to stream.
  virtual void write_body(FILE* fout) const;
  // Prints the profile in human-readable format to output stream.
  virtual void print(std::ostream& out) const;

 private:
  // Class identifier
  static const char* kClassID;

  // Return serialization class identity.
  virtual const char* class_id() const { return kClassID; }

  // Number of effective sequences in each alignment column.
  std::vector<float> neff_;
};  // CountProfile

}  // namespace cs

#endif  // SRC_COUNT_PROFILE_H_
