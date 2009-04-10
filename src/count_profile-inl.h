// Copyright 2009, Andreas Biegert

#ifndef SRC_COUNT_PROFILE_INL_H_
#define SRC_COUNT_PROFILE_INL_H_

#include "count_profile.h"

#include <string>
#include <vector>

namespace cs {

template<class Alphabet>
const char* CountProfile<Alphabet>::kClassID = "CountProfile";

template<class Alphabet>
inline CountProfile<Alphabet>::CountProfile(std::istream& in)
    : neff_(),
      has_counts_(false) {
  read(in);
}

template<class Alphabet>
inline CountProfile<Alphabet>::CountProfile(FILE* fin)
    : neff_(),
      has_counts_(false) {
  read(fin);
}

template<class Alphabet>
CountProfile<Alphabet>::CountProfile(const Sequence<Alphabet>& sequence)
    : Profile<Alphabet>(sequence.length()),
      neff_(sequence.length(), 1.0f),
      has_counts_(false) {
  for (int i = 0; i < num_cols(); ++i)
    data_[i][sequence[i]] = 1.0f;
}

template<class Alphabet>
CountProfile<Alphabet>::CountProfile(const Alignment<Alphabet>& alignment,
                                     bool position_specific_weights)
    : Profile<Alphabet>(alignment.num_match_cols()),
      neff_(alignment.num_match_cols()),
      has_counts_(false) {
  const int num_cols = alignment.num_match_cols();
  const int num_seqs = alignment.num_seqs();
  const int any   = Alphabet::instance().any();

  if (position_specific_weights) {
    matrix<float> w;  // position-specific sequence weights
    neff_ = position_specific_weights_and_diversity(alignment, w);
    for (int i = 0; i < num_cols; ++i)
      for (int k = 0; k < num_seqs; ++k)
        if (alignment[i][k] < any)
          data_[i][alignment[i][k]] += w[i][k];
  } else {
    std::vector<float> wg;  // global sequence weights
    neff_.assign(num_cols, global_weights_and_diversity(alignment, wg));
    for (int i = 0; i < num_cols; ++i)
      for (int k = 0; k < num_seqs; ++k)
        if (alignment[i][k] < any)
          data_[i][alignment[i][k]] += wg[k];
  }

  normalize(this);
}

template<class Alphabet>
inline CountProfile<Alphabet>::CountProfile(const CountProfile& other,
                                            int index,
                                            int length)
    : Profile<Alphabet>(other, index, length) {
  neff_.insert(neff_.begin(), other.neff_.begin() + index,
               other.neff_.begin() + index + length);
  has_counts_ = other.has_counts_;
}

template<class Alphabet>
void CountProfile<Alphabet>::readall(std::istream& in,
                                     std::vector< shared_ptr<CountProfile> >& v) {
  while (in.peek() && in.good()) {
    shared_ptr<CountProfile> p(new CountProfile(in));
    v.push_back(p);
  }
}

template<class Alphabet>
void CountProfile<Alphabet>::readall(FILE* fin,
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
void CountProfile<Alphabet>::convert_to_counts() {
  if (!has_counts_) {
    const bool islog = logspace();
    if (islog) transform_to_linspace();

    for (int i = 0; i < num_cols(); ++i)
      for (int a = 0; a < alphabet_size(); ++a)
        data_[i][a] *= neff_[i];
    has_counts_ = true;

    if (islog) transform_to_logspace();
  }
}

template<class Alphabet>
void CountProfile<Alphabet>::convert_to_frequencies() {
  if (has_counts_) {
    normalize(this);
    has_counts_ = false;
  }
}

template<class Alphabet>
void CountProfile<Alphabet>::read_header(std::istream& in) {
  Profile<Alphabet>::read_header(in);
  neff_.resize(num_cols());

  // Read has_counts
  std::string tmp;
  if (getline(in, tmp) && tmp.find("has_counts") != std::string::npos)
    has_counts_ = atoi(tmp.c_str() + 10) == 1;
}

template<class Alphabet>
void CountProfile<Alphabet>::read_header(FILE* fin) {
  Profile<Alphabet>::read_header(fin);
  neff_.resize(num_cols());

  // Read has_counts
  char buffer[kBufferSize];
  const char* ptr = buffer;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "has_counts")) {
    ptr = buffer;
    has_counts_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: profile does not contain 'has_counts' record!");
  }
}

template<class Alphabet>
void CountProfile<Alphabet>::read_body(std::istream& in) {
  std::string tmp;
  std::vector<std::string> tokens;
  int i = 0;
  getline(in, tmp);  // skip alphabet description line
  while (getline(in, tmp)) {
    if (tmp.empty()) continue;
    if (tmp.length() > 1 && tmp[0] == '/' && tmp[1] == '/') break;

    split(tmp, '\t', tokens);
    i = atoi(tokens[0].c_str()) - 1;
    for (int a = 0; a < alphabet_size(); ++a) {
      float log_p = tokens[a+1][0] == '*' ?
        std::numeric_limits<int>::max() : atoi(tokens[a+1].c_str());
      data_[i][a] =
        (logspace() ? -log_p / kLogScale : pow(2.0, -log_p / kLogScale));
    }
    neff_[i] = atof(tokens[alphabet_size()+1].c_str()) / kLogScale;
    tokens.clear();
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void CountProfile<Alphabet>::read_body(FILE* fin) {
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
        data_[i][a] = pow(2.0, static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
    }
    neff_[i] = static_cast<float>(strtoi(ptr)) / kLogScale;
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void CountProfile<Alphabet>::write_header(std::ostream& out) const {
  Profile<Alphabet>::write_header(out);
  out << "has_counts\t" << has_counts_ << std::endl;
}

template<class Alphabet>
void CountProfile<Alphabet>::write_body(std::ostream& out) const {
  out << "\t" << Alphabet::instance() << "\tNeff" << std::endl;
  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a) {
      float log_p = logspace() ? data_[i][a] : log2(data_[i][a]);
      if (-log_p == std::numeric_limits<float>::infinity())
        out << "\t*";
      else
        out << "\t" << -iround(log_p * kLogScale);
    }
    out << "\t" << iround(neff_[i] * kLogScale) << std::endl;
  }
  out << "//" << std::endl;
}

template<class Alphabet>
void CountProfile<Alphabet>::print(std::ostream& out) const {
  std::ios_base::fmtflags flags = out.flags();

  out << "\t" << Alphabet::instance() << "\tNeff" << std::endl;
  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a)
      out << '\t' << std::fixed << std::setprecision(4)
          << (logspace() ? pow(2.0, data_[i][a]) : data_[i][a]);
    out << '\t' << std::setprecision(2) << neff_[i] << std::endl;
  }

  out.flags(flags);
}

}  // namespace cs

#endif  // SRC_COUNT_PROFILE_INL_H_
