// Copyright 2009, Andreas Biegert

#ifndef SRC_COUNT_PROFILE_INL_H_
#define SRC_COUNT_PROFILE_INL_H_

#include "count_profile.h"

#include <vector>

#include "alignment-inl.h"
#include "exception.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* CountProfile<Alphabet>::kClassID = "CountProfile";

template<class Alphabet>
inline CountProfile<Alphabet>::CountProfile(FILE* fin)
    : Profile<Alphabet>(),
      neff_() {
  Read(fin);
}

template<class Alphabet>
CountProfile<Alphabet>::CountProfile(const Sequence<Alphabet>& sequence)
    : Profile<Alphabet>(sequence.length()),
      neff_(sequence.length(), 1.0f) {
  for (int i = 0; i < num_cols(); ++i)
    data_[i][sequence[i]] = 1.0f;
}

template<class Alphabet>
CountProfile<Alphabet>::CountProfile(const Alignment<Alphabet>& alignment,
                                     bool position_specific_weights)
    : Profile<Alphabet>(alignment.num_match_cols()),
      neff_(alignment.num_match_cols()) {
  const int num_cols = alignment.num_match_cols();
  const int num_seqs = alignment.num_seqs();
  const int any      = Alphabet::instance().any();

  if (position_specific_weights) {
    matrix<float> w;  // position-specific sequence weights
    neff_ = PositionSpecificWeightsAndDiversity(alignment, w);
    for (int i = 0; i < num_cols; ++i) {
      for (int k = 0; k < num_seqs; ++k)
        if (alignment[i][k] < any)
          data_[i][alignment[i][k]] += w[i][k];
      data_[i][any] = 1.0f;
    }
  } else {
    std::vector<float> wg;  // global sequence weights
    neff_.assign(num_cols, GlobalWeightsAndDiversity(alignment, wg));
    for (int i = 0; i < num_cols; ++i) {
      for (int k = 0; k < num_seqs; ++k)
        if (alignment[i][k] < any)
          data_[i][alignment[i][k]] += wg[k];
      data_[i][any] = 1.0f;
    }
  }

  Normalize(this);
}

template<class Alphabet>
inline CountProfile<Alphabet>::CountProfile(const CountProfile& other,
                                            int index,
                                            int length)
    : Profile<Alphabet>(other, index, length) {
  neff_.insert(neff_.begin(), other.neff_.begin() + index,
               other.neff_.begin() + index + length);
}

template<class Alphabet>
void CountProfile<Alphabet>::ReadAll(
    FILE* fin,
    std::vector< shared_ptr<CountProfile> >* v) {
  while (!feof(fin)) {
    shared_ptr<CountProfile> p(new CountProfile(fin));
    v->push_back(p);

    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
}

template<class Alphabet>
void CountProfile<Alphabet>::ReadHeader(FILE* fin) {
  Profile<Alphabet>::ReadHeader(fin);
  neff_.resize(num_cols());
}

template<class Alphabet>
void CountProfile<Alphabet>::ReadBody(FILE* fin) {
  const int alph_size = alphabet_size();
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int i = 0;

  fgetline(buffer, kBufferSize, fin);  // skip alphabet description line
  while (fgetline(buffer, kBufferSize, fin)
         && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    i = strtoi(ptr) - 1;
    for (int a = 0; a < alph_size; ++a) {
      if (logspace())
        data_[i][a] = static_cast<float>(-strtoi_ast(ptr)) / kLogScale;
      else
        data_[i][a] = fast_pow2(static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
    }
    data_[i][alph_size] = logspace() ? 0.0f : 1.0f;
    neff_[i] = static_cast<float>(strtoi(ptr)) / kLogScale;
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void CountProfile<Alphabet>::WriteBody(FILE* fout) const {
  fputs("PROF\t", fout);
  Alphabet::instance().Write(fout);
  fputs("\tNEFF\n", fout);

  for (int i = 0; i < num_cols(); ++i) {
    fprintf(fout, "%i", i+1);
    for (int a = 0; a < alphabet_size(); ++a) {
      float log_p = logspace() ? data_[i][a] : fast_log2(data_[i][a]);
      if (log_p == -INFINITY)
        fputs("\t*", fout);
      else
        fprintf(fout, "\t%i", -iround(log_p * kLogScale));
    }
    fprintf(fout, "\t%i\n", iround(neff_[i] * kLogScale));
  }
  fputs("//\n", fout);
}

template<class Alphabet>
void CountProfile<Alphabet>::Print(std::ostream& out) const {
  out << "\t" << Alphabet::instance() << "\tNeff" << std::endl;

  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a)
      out << strprintf("\t%6.4f",
                       logspace() ? fast_pow2(data_[i][a]) : data_[i][a]);
    out << strprintf("\t%-5.2f\n", neff_[i]);
  }
}

}  // namespace cs

#endif  // SRC_COUNT_PROFILE_INL_H_
