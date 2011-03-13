// Copyright 2009, Andreas Biegert

#ifndef CS_TRAINING_SEQUENCE_H_
#define CS_TRAINING_SEQUENCE_H_

#include "profile_column.h"
#include "sequence-inl.h"

namespace cs {

template<class Abc>
struct TrainingSequence {
  TrainingSequence(size_t wlen) : x(wlen) {}

  TrainingSequence(const Sequence<Abc>& seq, const ProfileColumn<Abc>& col)
    : x(seq), y(col) {}

  TrainingSequence(FILE* fin) {
    if (!StreamStartsWith(fin, "TrainingSequence"))
      throw Exception("Stream does not start with class id 'TrainingSequence'!");
    char buffer[KB];
    fgetline(buffer, KB, fin);  // skip alphabet description line
    fgetline(buffer, KB, fin);
    const char* ptr = buffer;
    for (size_t a = 0; a < Abc::kSizeAny; ++a)
      y[a] = exp(static_cast<double>(-strastoi(ptr)) / kScale);
    x.Read(fin);
  }

  void Write(FILE* fout) const {
    fputs("TrainingSequence\n", fout);
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

  Sequence<Abc> x;       // input sequence window representing context
  ProfileColumn<Abc> y;  // pseudocounts at central column to be predicted
};  // class TrainingSequence

// Prints training sequence for debugging
template<class Abc>
std::ostream& operator<< (std::ostream& out, const TrainingSequence<Abc>& tseq) {
  out << "TrainingSequence" << std::endl << tseq.y << tseq.x;
  return out;
}

}  // namespace cs

#endif  // CS_TRAINING_SEQUENCE_H_
