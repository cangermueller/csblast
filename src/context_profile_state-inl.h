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

#ifndef CS_CONTEXT_PROFILE_STATE_INL_H_
#define CS_CONTEXT_PROFILE_STATE_INL_H_

#include "context_profile_state.h"

#include <stdio.h>
#include <string.h>

#include "context_profile-inl.h"
#include "profile-inl.h"
#include "transition.h"

namespace cs {

template<class Alphabet>
const char* ContextProfileState<Alphabet>::kClassID = "ContextProfileState";

template<class Alphabet>
ContextProfileState<Alphabet>::ContextProfileState(FILE* fin)
    : ContextProfile<Alphabet>() {
  Read(fin);
}

template<class Alphabet>
ContextProfileState<Alphabet>::ContextProfileState(
    int index,
    int num_states,
    const Profile<Alphabet>& profile)
    : ContextProfile<Alphabet>(index, profile),
      in_transitions_(num_states),
      out_transitions_(num_states) {}

template<class Alphabet>
void ContextProfileState<Alphabet>::ClearTransitions() {
  in_transitions_.clear();
  out_transitions_.clear();
}

template<class Alphabet>
void ContextProfileState<Alphabet>::ReadHeader(FILE* fin) {
  ContextProfile<Alphabet>::ReadHeader(fin);

  // Read HMM size
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int num_states = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'NSTATES' record!");
  }

  in_transitions_.resize(num_states);
  out_transitions_.resize(num_states);
}

template<class Alphabet>
void ContextProfileState<Alphabet>::WriteHeader(FILE* fout) const {
  ContextProfile<Alphabet>::WriteHeader(fout);
  fprintf(fout, "NSTATES\t%i\n", static_cast<int>(in_transitions_.size()));
}

}  // namespace cs

#endif  // CS_CONTEXT_PROFILE_STATE_INL_H_
