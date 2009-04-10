// Copyright 2009, Andreas Biegert

#ifndef SRC_SEQUENCE_INL_H_
#define SRC_SEQUENCE_INL_H_

#include "sequence.h"

#include <string>
#include <vector>

namespace cs {

template<class Alphabet>
inline Sequence<Alphabet>::Sequence(int length)
    : seq_(length) { }

template<class Alphabet>
inline Sequence<Alphabet>::Sequence(std::istream& in) {
  read(in);
}

template<class Alphabet>
inline Sequence<Alphabet>::Sequence(FILE* in) {
  read(in);
}

template<class Alphabet>
Sequence<Alphabet>::Sequence(const std::string& header,
                             const std::string& sequence) {
  init(header, sequence);
}

template<class Alphabet>
inline void Sequence<Alphabet>::readall(std::istream& in,
                                        std::vector< shared_ptr<Sequence> >& v) {
  while (in.peek() && in.good()) {
    shared_ptr<Sequence> p(new Sequence(in));
    v.push_back(p);
  }
}

template<class Alphabet>
inline void Sequence<Alphabet>::readall(FILE* fin,
                                        std::vector< shared_ptr<Sequence> >* v) {
  while (!feof(fin)) {
    shared_ptr<Sequence> p(new Sequence(fin));
    v->push_back(p);
  }
}

template<class Alphabet>
void Sequence<Alphabet>::init(std::string header, std::string sequence) {
  // Init header
  std::string().swap(header_);
  header_.append(header.begin() + (header[0] == '>' ? 1 : 0), header.end());

  // Strip whitespace and newlines from sequence.
  sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace),
                 sequence.end());
  // Validate each character and convert to integer representation
  const Alphabet& alphabet = Alphabet::instance();
  const int seqlen = sequence.length();
  seq_.resize(seqlen);
  for (int i = 0; i < seqlen; ++i) {
    char c = sequence[i];
    if (alphabet.valid(c)) {
      seq_[i] = alphabet.ctoi(c);
    } else {
      throw Exception("Invalid character %c at position %i of sequence '%s'",
                      c, i, header_.c_str());
    }
  }
}

template<class Alphabet>
void Sequence<Alphabet>::read(std::istream& in) {
  LOG(DEBUG1) << "Reading sequence from stream ...";
  std::string tmp;
  std::string header;
  std::string sequence;

  // Read header
  if (getline(in, tmp)) {
    if (tmp.empty() ||  tmp[0] != '>')
      throw Exception("Sequence header does not start with '>' character!");
    header.append(tmp.begin() + 1, tmp.end());
  } else {
    throw Exception("Failed to read sequence from stream!");
  }
  // Read sequence
  while (in.peek() != '>' && getline(in, tmp))
    sequence.append(tmp.begin(), tmp.end());

  init(header, sequence);
  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void Sequence<Alphabet>::read(FILE* fin) {
  LOG(DEBUG1) << "Reading sequence from stream ...";

  char buffer[kBufferSize];
  int c = '\0';
  std::string header;
  std::string sequence;

  // Read header
  while (fgetline(buffer, kBufferSize, fin)) {
    if (!strscn(buffer)) continue;
    if (buffer[0] == '>') {
      header.append(buffer + 1);
      break;
    } else {
      throw Exception("Sequence header does not start with '>'!");
    }
  }
  // Read sequence
  while (fgetline(buffer, kBufferSize, fin)) {
    if (strscn(buffer))
      sequence.append(buffer);

    c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
    if (static_cast<char>(c) == '>') break;
  }
  init(header, sequence);

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void Sequence<Alphabet>::write(std::ostream& out, int width) const {
  out << '>' << header_ << std::endl;
  for (int i = 0; i < length(); ++i) {
    out << chr(i);
    if ((i+1) % width == 0) out << std::endl;
  }
  if (length() % width != 0) out << std::endl;
}

}  // namespace cs

#endif  // SRC_SEQUENCE_INL_H_
