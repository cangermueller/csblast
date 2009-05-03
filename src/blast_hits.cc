// Copyright 2009, Andreas Biegert

#include "blast_hits.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>

#include <string>
#include <vector>

#include "globals.h"
#include "exception.h"
#include "log.h"
#include "utils-inl.h"

using std::string;
using std::vector;

namespace cs {

BlastHits::BlastHits() {}

BlastHits::BlastHits(FILE* fin) {
  Read(fin);
}

void BlastHits::Read(FILE* fin) {
  char buffer[KB];
  const char* ptr;
  hits_.clear();

  // Advance to hitlist and parse query length on the way
  while (fgetline(buffer, KB, fin)) {
    if (strstr(buffer, "Sequences producing significant alignments")) break;
    if (strstr(buffer, "letters)")) {
      ptr = buffer;
      query_length_ = strtoi(ptr);
    }
  }
  fgetline(buffer, KB, fin);  // skip empty line

  // Parse hitlist
  while (fgetline(buffer, KB, fin)) {
    if (!strscn(buffer)) break;  // reached end of hitlist
    ptr = buffer + strlen(buffer) - 1;

    while (isspace(*ptr)) --ptr;  // find end of e-value
    if (!isdigit(*ptr)) break;    // broken hit without bit-score and evalue

    Hit h;
    while (isgraph(*ptr)) --ptr;  // find start of e-value
    h.evalue = atof(ptr + 1);
    while (isspace(*ptr)) --ptr;  // find end of bit-score
    while (isgraph(*ptr)) --ptr;  // find start of bit-score
    h.bit_score = atof(ptr + 1);
    while (isspace(*ptr)) --ptr;  // find end of defline
    h.definition = string(buffer, ptr - buffer + 1);
    h.oid = hits_.size() + 1;
    hits_.push_back(h);
  }

  // Parse alignments
  int i = -1;  // hit index
  while (fgetline(buffer, KB, fin)) {
    if (!strscn(buffer)) continue;

    if (strstr(buffer, "Database:")) {
      break;

    } else if (buffer[0] == '>') {
      ++i;
      if (i >= static_cast<int>(hits_.size())) break;

    } else if (strstr(buffer, "Score =")) {
      Hsp hsp;
      ptr = strchr(buffer, '=') + 1;
      while (isspace(*ptr)) ++ptr;
      hsp.bit_score = atof(ptr);
      ptr = strchr(ptr, '=') + 1;
      while (isspace(*ptr)) ++ptr;
      hsp.evalue = atof(ptr);
      hits_[i].hsps.push_back(hsp);

    } else if (strstr(buffer, "Query:")) {
      ptr = buffer;
      int query_start = strtoi(ptr);
      if (hits_[i].hsps.back().query_start == 0)
        hits_[i].hsps.back().query_start = query_start;
      while (isspace(*ptr)) ++ptr;
      while (isgraph(*ptr)) hits_[i].hsps.back().query_seq.push_back(*ptr++);
      hits_[i].hsps.back().query_end = strtoi(ptr);
      hits_[i].hsps.back().length = hits_[i].hsps.back().query_seq.size();

    } else if (strstr(buffer, "Sbjct:")) {
      ptr = buffer;
      int subject_start = strtoi(ptr);
      if (hits_[i].hsps.back().subject_start == 0)
        hits_[i].hsps.back().subject_start = subject_start;
      while (isspace(*ptr)) ++ptr;
      while (isgraph(*ptr)) hits_[i].hsps.back().subject_seq.push_back(*ptr++);
      hits_[i].hsps.back().subject_end = strtoi(ptr);
    }
  }
}

void BlastHits::Filter(double evalue_threshold) {
  // Delete hits below E-value threshold
  for (HitIter hit = begin(); hit != end(); ++hit) {
    if (hit->evalue > evalue_threshold) {
      hits_.erase(hit, end());
      break;
    }
  }

  // Delete HSPs in remaining hits with E-value below threshold
  for (HitIter hit = begin(); hit != end(); ++hit) {
    for (HspIter hsp = hit->hsps.begin(); hsp != hit->hsps.end(); ++hsp) {
      if (hsp->evalue > evalue_threshold) {
        hit->hsps.erase(hsp, hit->hsps.end());
        break;
      }
    }
  }
}

}  // namespace cs
