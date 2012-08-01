/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_PROFILE_LIBRARY_INL_H_
#define CS_PROFILE_LIBRARY_INL_H_

#include "profile_library.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils.h"

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
  Read(fin);
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
  profile_init.Init(*this);
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::clear() {
  profiles_.clear();
  profiles_.reserve(num_profiles());
}

template<class Alphabet>
inline int ProfileLibrary<Alphabet>::AddProfile(
    const Profile<Alphabet>& profile) {
  assert_eq(num_cols(), profile.num_cols());
  if (full()) return -1;

  shared_ptr< ContextProfile<Alphabet> >
    profile_ptr(new ContextProfile<Alphabet>(profiles_.size(), profile));
  profile_ptr->set_prior(1.0f / num_profiles());

  profiles_.push_back(profile_ptr);
  return profiles_.size() - 1;
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::TransformToLogSpace() {
  if (!logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->TransformToLogSpace();
    logspace_ = true;
  }
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::TransformToLinSpace() {
  if (logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->TransformToLinSpace();
    logspace_ = false;
  }
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::Read(FILE* fin) {
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
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NPROF")) {
    ptr = buffer;
    num_profiles_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'NPROF' record!");
  }
  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NCOLS")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'NCOLS' record!");
  }
  // Read number of iterations
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ITERS")) {
    ptr = buffer;
    iterations_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'ITERS' record!");
  }
  // Read states logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "LOG")) {
    ptr = buffer;
    logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: HMM does not contain 'LOG' record!");
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
void ProfileLibrary<Alphabet>::Write(FILE* fout) const {
  // Write header
  fputs("ProfileLibrary\n", fout);
  fprintf(fout, "NPROF\t%i\n", num_profiles());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());
  fprintf(fout, "LOG\t%i\n", logspace() ? 1 : 0);

  // Serialize profiles
  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    (*pi)->Write(fout);
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::Print(std::ostream& out) const {
  out << "ProfileLibrary" << std::endl;
  out << "Total number of profiles: " << num_profiles() << std::endl;
  out << "Context profile columns:  " << num_cols() << std::endl;
  out << "Clustering iterations:    " << iterations() << std::endl;

  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    out << **pi;
}


template<class Alphabet>
void SamplingProfileInitializer<Alphabet>::Init(
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
    if (pc_) pc_->AddPseudocountsToProfile(ConstantAdmixture(pc_admixture_), &p);
    lib.AddProfile(p);
  }
  if (!lib.full())
    throw Exception("Could not fully initialize all %i library profiles. "
                    "Maybe too few training profiles provided?",
                    lib.num_profiles());

  LOG(DEBUG) << "Profile library after profile initialization:";
  LOG(DEBUG) << lib;
}

}  // namespace cs

#endif  // CS_PROFILE_LIBRARY_INL_H_
