// Copyright 2009, Andreas Biegert

#ifndef CS_TRAINING_PROFILE_H_
#define CS_TRAINING_PROFILE_H_

#include "count_profile-inl.h"
#include "profile_column.h"

namespace cs {

template<class Abc>
struct TrainingProfile {
  TrainingProfile(size_t wlen) : x(wlen) {}

  TrainingProfile(const Sequence<Abc>& seq, const ProfileColumn<Abc>& col)
    : x(seq), y(col) {}

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

  void Write(FILE* fout) const {
    fputs("TrainingProfile\n", fout);
    x.Write(fout);
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
    fputs("//\n", fout);
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
