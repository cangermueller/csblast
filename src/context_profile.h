// Copyright 2009, Andreas Biegert

#ifndef SRC_CONTEXT_PROFILE_H_
#define SRC_CONTEXT_PROFILE_H_

#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "abstract_state.h"
#include "alignment.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs {

// A container class for profiles derived from alignments.
template<class Alphabet>
class ContextProfile : public Profile<Alphabet> {
 public:
  // Needed to access names in templatized Profile base class
  using Profile<Alphabet>::num_cols;
  using Profile<Alphabet>::alphabet_size;
  using Profile<Alphabet>::logspace;
  using Profile<Alphabet>::Read;

  // Constructs a dummy context profile.
  ContextProfile();
  // Constructs a context profile with num_cols columns initialized to zero.
  ContextProfile(int index, int num_cols);
  // Constructs profile from serialized profile read from input stream.
  explicit ContextProfile(FILE* fin);
  // Constructs a context profile from simple profile and checks if length is
  // valid.
  ContextProfile(int index, const Profile<Alphabet>& profile);

  virtual ~ContextProfile() {}

  // Returns the index of this context profile.
  int index() const { return index_; }
  // Sets the index of this context profile.
  void set_index(int i) { index_ = i; }
  // Returns the prior probability of this context profile.
  double prior() const { return prior_; }
  // Sets the prior probability of this context profile.
  void set_prior(double p) { prior_ = p; }
  // Returns index of central profile column.
  int center() const { return (num_cols() - 1) / 2; }

 protected:
  // Needed to access names in templatized Profile base class
  using Profile<Alphabet>::data_;
  using Profile<Alphabet>::kBufferSize;

  // Reads and initializes serialized scalar data members from stream.
  virtual void ReadHeader(FILE* fin);
  // Writes serialized scalar data members to stream.
  virtual void WriteHeader(FILE* fout) const;
  // Prints the profile in human-readable format to output stream.
  virtual void Print(std::ostream& out) const;

  // Index of context-profile.
  int index_;
  // Prior probability of context-profile.
  double prior_;

 private:
  // Class identifier
  static const char* kClassID;

  // Return serialization class identity.
  virtual const char* class_id() const { return kClassID; }
  // Checks if profile has odd number of columns.
  void check();
};  // ContextProfile

}  // namespace cs

#endif  // SRC_CONTEXT_PROFILE_H_
