/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

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

#ifndef CS_TRAINING_PROFILE_H_
#define CS_TRAINING_PROFILE_H_

#include "count_profile-inl.h"
#include "profile_column.h"
#include "training_sequence.h"

namespace cs {

template<class Abc>
struct TrainingProfile {
  TrainingProfile(size_t wlen) : x(wlen) {}

  TrainingProfile(const Sequence<Abc>& seq, const ProfileColumn<Abc>& col)
    : x(seq), y(col) {}

  TrainingProfile(const CountProfile<Abc>& cp, const ProfileColumn<Abc>& col)
    : x(cp), y(col) {}

  TrainingProfile(const TrainingSequence<Abc>& tseq)
    : x(CountProfile<Abc>(tseq.x)), y(tseq.y) {}

  TrainingProfile(FILE* fin) {
    if (!StreamStartsWith(fin, "TrainingProfile"))
      throw Exception("Stream does not start with class id 'TrainingProfile'!");
    char buffer[KB];
    fgetline(buffer, KB, fin);  // skip alphabet description line
    fgetline(buffer, KB, fin);
    const char* ptr = buffer;
    for (size_t a = 0; a < Abc::kSize; ++a)
      y[a] = exp(static_cast<double>(-strastoi(ptr)) / kScale);
    x.Read(fin);
  }

  static bool IsTrainingProfile(FILE* fin) {
    return StreamStartsWith(fin, "TrainingProfile");
  }

  void Write(FILE* fout) const {
    fputs("TrainingProfile\n", fout);
    for (size_t a = 0; a < Abc::kSizeAny; ++a)
      fprintf(fout, "\t%c", Abc::kIntToChar[a]);
    fputs("\nPC", fout);
    for (size_t a = 0; a < Abc::kSizeAny; ++a) {
      double log_val = log(y[a]);
      if (log_val == -INFINITY) fputs("\t*", fout);
      else fprintf(fout, "\t%d", -iround(log_val * kScale));
    }
    fputs("\n", fout);
    x.Write(fout);
  }

  CountProfile<Abc> x;   // input count profile window representing context
  ProfileColumn<Abc> y;  // pseudocounts at central column to be predicted
};  // class TrainingProfile

// Prints training profile for debugging
template<class Abc>
std::ostream& operator<< (std::ostream& out, const TrainingProfile<Abc>& tprof) {
  out << "TrainingProfile" << std::endl << tprof.y << tprof.x;
  return out;
}

}  // namespace cs

#endif  // CS_TRAINING_PROFILE_H_
