// Copyright 2009, Andreas Biegert

#ifndef SRC_PROFILE_LIBRARY_INL_H_
#define SRC_PROFILE_LIBRARY_INL_H_

#include "profile_library.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>

#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(int num_profiles, int num_cols)
    : num_profiles_(num_profiles),
      num_cols_(num_cols),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  profiles_.reserve(num_profiles);
}

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(FILE* fin)
    : num_profiles_(0),
      num_cols_(0),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  read(fin);
}

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(
    int num_profiles,
    int num_cols,
    const ProfileInitializer<Alphabet>& profile_init)
    : num_profiles_(num_profiles),
      num_cols_(num_cols),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  profiles_.reserve(num_profiles);
  profile_init.init(*this);
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::clear() {
  profiles_.clear();
  profiles_.reserve(num_profiles());
}

template<class Alphabet>
inline int ProfileLibrary<Alphabet>::add_profile(
    const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("Profile library contains already %i profiles!",
                    num_profiles());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  shared_ptr< ContextProfile<Alphabet> >
    profile_ptr(new ContextProfile<Alphabet>(profiles_.size(), profile));
  profile_ptr->set_prior(1.0f / num_profiles());

  profiles_.push_back(profile_ptr);
  return profiles_.size() - 1;
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::transform_to_logspace() {
  if (!logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->transform_to_logspace();
    logspace_ = true;
  }
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::transform_to_linspace() {
  if (logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->transform_to_linspace();
    logspace_ = false;
  }
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::read(FILE* fin) {
  LOG(DEBUG1) << "Reading profile library from stream ...";

  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Check if stream actually contains a serialized HMM
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, "ProfileLibrary"))
    throw Exception("Profile library does not start with 'ProfileLibrary' "
                    "keyword!");

  // Read number of profiles
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "num_profiles")) {
    ptr = buffer;
    num_profiles_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'num_profiles' record!");
  }
  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "num_cols")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'num_cols' record!");
  }
  // Read number of iterations
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "iterations")) {
    ptr = buffer;
    iterations_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'iterations' record!");
  }
  // Read states logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "logspace")) {
    ptr = buffer;
    logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: HMM does not contain 'logspace' record!");
  }

  // Read context profiles
  profiles_.reserve(num_profiles());
  while (!full() && !feof(fin)) {
    shared_ptr< ContextProfile<Alphabet> > p(new ContextProfile<Alphabet>(fin));
    profiles_.push_back(p);
  }
  if (!full())
    throw Exception("Profile library has %i profiles but should have %i!",
                    profiles_.size(), num_profiles());

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::write(FILE* fout) const {
  // Write header
  fputs("ProfileLibrary\n", fout);
  fprintf(fout, "num_profiles\t%i\n", num_profiles());
  fprintf(fout, "num_cols\t%i\n", num_cols());
  fprintf(fout, "iterations\t%i\n", iterations());
  fprintf(fout, "logspace\t%i\n", logspace() ? 1 : 0);

  // Serialize profiles
  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    (*pi)->write(fout);
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::print(std::ostream& out) const {
  out << "ProfileLibrary" << std::endl;
  out << "Total number of profiles: " << num_profiles() << std::endl;
  out << "Context profile columns:  " << num_cols() << std::endl;
  out << "Clustering iterations:    " << iterations() << std::endl;

  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    out << **pi;
}


template<class Alphabet>
void SamplingProfileInitializer<Alphabet>::init(
    ProfileLibrary<Alphabet>& lib) const {
  LOG(DEBUG) << "Initializing profile library with " << lib.num_profiles()
             << " profile windows randomly sampled from "
             << profiles_.size() << " training profiles ...";

  // Iterate over randomly shuffled profiles and add them to the library until
  // the library is full.
  for (profile_iterator pi = profiles_.begin();
       pi != profiles_.end() && !lib.full(); ++pi) {
    if (lib.num_cols() != (*pi)->num_cols())
      throw Exception("Library has num_cols=%i but training profile has %i!",
                      lib.num_cols(), (*pi)->num_cols());

    CountProfile<Alphabet> p(**pi);
    p.convert_to_frequencies();
    if (pc_) pc_->add_to_profile(ConstantAdmixture(pc_admixture_), &p);
    lib.add_profile(p);
  }
  if (!lib.full())
    throw Exception("Could not fully initialize all %i library profiles. "
                    "Maybe too few training profiles provided?",
                    lib.num_profiles());

  LOG(DEBUG) << "Profile library after profile initialization:";
  LOG(DEBUG) << lib;
}

}  // namespace cs

#endif  // SRC_PROFILE_LIBRARY_INL_H_
