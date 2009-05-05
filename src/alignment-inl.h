// Copyright 2009, Andreas Biegert

#ifndef SRC_ALIGNMENT_INL_H_
#define SRC_ALIGNMENT_INL_H_

#include "alignment.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "blast_hits.h"
#include "exception.h"
#include "log.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
Alignment<Alphabet>::Alignment(FILE* fin, Format format) {
  read(fin, format);
}

template<class Alphabet>
Alignment<Alphabet>::Alignment(const Sequence<Alphabet>& seq) {
  std::vector<std::string> headers;
  std::vector<std::string> seqs;
  headers.push_back(seq.header());
  seqs.push_back(seq.ToString());
  Init(headers, seqs);
}

template<class Alphabet>
Alignment<Alphabet>::Alignment(const BlastHits& hits, bool best) {
  std::vector<std::string> headers;
  std::vector<std::string> seqs;

  typedef typename BlastHits::ConstHitIter HitIter;
  typedef typename BlastHits::ConstHspIter HspIter;
  for (HitIter hit = hits.begin(); hit != hits.end(); ++hit) {
    for (HspIter hsp = hit->hsps.begin(); hsp != hit->hsps.end(); ++hsp) {
      // Construct query anchored alignment string
      std::string seq(hsp->query_start - 1, '-');
      for (int i =  0; i < hsp->length; ++i)
        if (hsp->query_seq[i] != '-')
          seq += hsp->subject_seq[i];
      seq.append(hits.query_length() - seq.length(), '-');

      headers.push_back(hit->definition);
      seqs.push_back(seq);

      if (best) break;
    }
  }
  Init(headers, seqs);
}

template<class Alphabet>
void Alignment<Alphabet>::readall(FILE* fin,
                                  Format format,
                                  std::vector< shared_ptr<Alignment> >* v) {
  while (!feof(fin)) {
    shared_ptr<Alignment> p(new Alignment(fin, format));
    v->push_back(p);

    int c = fgetc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
}

template<class Alphabet>
void Alignment<Alphabet>::Init(const std::vector<std::string>& headers,
                               const std::vector<std::string>& seqs) {
  if (seqs.empty())
    throw Exception("Bad alignment: initialization with empty sequence vector!");
  if (headers.size() != seqs.size())
    throw Exception("Bad alignment: unequal number of headers and sequences!");

  const int num_seqs = seqs.size();
  const int num_cols = seqs[0].length();
  for (int k = 1; k < num_seqs; ++k)
    if (static_cast<int>(seqs[k].length()) != num_cols)
      throw Exception("Alignment sequence %i has length %i but should have %i!",
                      k+1, seqs[k].length(), num_cols);

  // Validate characters and convert to integer representation
  const Alphabet& alphabet = Alphabet::instance();
  headers_.resize(num_seqs);
  resize(seqs.size(), seqs[0].length());
  for (int k = 0; k < num_seqs; ++k) {
    headers_[k] = headers[k];
    for (int i = 0; i < num_cols; ++i) {
      const char c = seqs[k][i];
      column_indexes_[i] = i;
      match_column_[i] = match_chr(c);
      if (alphabet.valid(c, true))
        seqs_[i][k] = alphabet.ctoi(c);
      else
        throw Exception("Invalid character %c at position %i of sequence '%s'",
                        c, i, headers_[k].c_str());
    }
  }

  // Replace gap with endgap for all gaps at either end of a sequence
  for (int k = 0; k < num_seqs; ++k) {
    for (int i = 0; i < num_cols && seqs_[i][k] == alphabet.gap(); ++i)
      seqs_[i][k] = alphabet.endgap();
    for (int i = num_cols - 1; i >= 0 && seqs_[i][k] == alphabet.gap(); --i)
      seqs_[i][k] = alphabet.endgap();
  }

  // Initialize index array for match columns
  set_match_indexes();
}

template<class Alphabet>
void Alignment<Alphabet>::set_match_indexes() {
  const int match_cols = std::count(&match_column_[0],
                                    &match_column_[0] + num_cols(), true);
  match_indexes.resize(match_cols);
  match_indexes = column_indexes_[match_column_];
}

template<class Alphabet>
void Alignment<Alphabet>::read(FILE* fin, Format format) {
  LOG(DEBUG4) << "Reading alignment from stream ...";

  std::vector<std::string> headers;
  std::vector<std::string> seqs;
  switch (format) {
    case FASTA:
      read_fasta(fin, &headers, &seqs);
      break;
    case A2M:
      read_a2m(fin, &headers, &seqs);
      break;
    case A3M:
      read_a3m(fin, &headers, &seqs);
      break;
    default:
      throw Exception("Unsupported alignment input format %i!", format);
  }
  Init(headers, seqs);

  LOG(DEBUG4) << *this;
}

template<class Alphabet>
void Alignment<Alphabet>::read_fasta_flavors(FILE* fin,
                                             std::vector<std::string>* headers,
                                             std::vector<std::string>* seqs) {
  headers->clear();
  seqs->clear();

  char buffer[kBufferSize];
  int c = '\0';
  while (!feof(fin) && static_cast<char>(c) != '#') {
    // Read header
    while (fgetline(buffer, kBufferSize, fin)) {
      if (!strscn(buffer)) continue;
      if (buffer[0] == '>') {
        headers->push_back(std::string(buffer + 1));
        break;
      } else {
        throw Exception("Header of sequence %i does not start with '>'!",
                        headers->size() + 1);
      }
    }

    // Read sequence
    seqs->push_back("");
    while (fgetline(buffer, kBufferSize, fin)) {
      if (!strscn(buffer)) continue;
      seqs->back().append(buffer);

      c = getc(fin);
      if (c == EOF || static_cast<char>(c) == '#') break;
      ungetc(c, fin);
      if (static_cast<char>(c) == '>') break;
    }
    // Remove whitespace
    seqs->back().erase(remove_if(seqs->back().begin(), seqs->back().end(), isspace),
                       seqs->back().end());

    LOG(DEBUG2) << headers->back();
  }
  if (headers->empty())
    throw Exception("Bad alignment: no alignment data found in stream!");
}


template<class Alphabet>
void Alignment<Alphabet>::read_fasta(FILE* fin,
                                     std::vector<std::string>* headers,
                                     std::vector<std::string>* seqs) {
  read_fasta_flavors(fin, headers, seqs);
  // convert all characters to match characters
  for (std::vector<std::string>::iterator it = seqs->begin();
       it != seqs->end(); ++it)
    transform(it->begin(), it->end(), it->begin(), to_match_chr);
}

template<class Alphabet>
void Alignment<Alphabet>::read_a2m(FILE* fin,
                                   std::vector<std::string>* headers,
                                   std::vector<std::string>* seqs) {
  read_fasta_flavors(fin, headers, seqs);
}

template<class Alphabet>
void Alignment<Alphabet>::read_a3m(FILE* fin,
                                   std::vector<std::string>* headers,
                                   std::vector<std::string>* sequences) {
  read_fasta_flavors(fin, headers, sequences);

  // Check number of match states
  std::vector<std::string>& seqs    = *sequences;
  const int num_seqs = seqs.size();
  const int num_match_cols =
    count_if(seqs[0].begin(), seqs[0].end(), match_chr);
  for (int k = 1; k < num_seqs; ++k) {
    const int num_match_cols_k =
      count_if(seqs[k].begin(), seqs[k].end(), match_chr);
    if (num_match_cols_k != num_match_cols)
      throw Exception("Sequence %i has %i match columns but should have %i!",
                      k, num_match_cols_k, num_match_cols);
    if (count(seqs[k].begin(), seqs[k].end(), '.') > 0)
      throw Exception("Sequence %i in A3M alignment contains gaps!", k);
  }

  // Insert gaps into A3M alignment
  std::vector<std::string> seqs_a2m(seqs.size(), "");
  matrix<std::string> inserts(seqs.size(), num_match_cols, "");
  std::vector<int> max_insert_len(num_match_cols, 0);

  // Move inserts before first match state into seqs_a2m and keep track of
  // longest first insert
  int max_first_insert_len = 0;
  for (int k = 0; k < num_seqs; ++k) {
    std::string::iterator i =
      find_if(seqs[k].begin(), seqs[k].end(), match_chr);
    if (i != seqs[k].end()) {
      seqs_a2m[k].append(seqs[k].begin(), i);
      seqs[k].erase(seqs[k].begin(), i);
    }
    max_first_insert_len = std::max(max_first_insert_len,
                                    static_cast<int>(i - seqs[k].begin()));
  }

  // Extract all inserts and keep track of longest insert after each match
  // column
  for (int k = 0; k < num_seqs; ++k) {
    int i = -1;
    std::string::iterator insert_end = seqs[k].begin();
    std::string::iterator insert_start =
      find_if(insert_end, seqs[k].end(), insert_chr);
    while (insert_start != seqs[k].end() && insert_end != seqs[k].end()) {
      i += insert_start - insert_end;
      insert_end = find_if(insert_start, seqs[k].end(), match_chr);
      inserts[k][i] = std::string(insert_start, insert_end);
      max_insert_len[i] = std::max(static_cast<int>(inserts[k][i].length()),
                                   max_insert_len[i]);
      insert_start = find_if(insert_end, seqs[k].end(), insert_chr);
    }
  }

  // Build new A2M alignment
  for (int k = 0; k < num_seqs; ++k) {
    seqs_a2m[k].append(max_first_insert_len - seqs_a2m[k].length(), '.');
    int i = 0;
    std::string::iterator match = seqs[k].begin();
    while (match != seqs[k].end()) {
      seqs_a2m[k].append(1, *match);
      if (max_insert_len[i] > 0) {
        seqs_a2m[k].append(inserts[k][i]);
        seqs_a2m[k].append(max_insert_len[i] - inserts[k][i].length(), '.');
      }
      match = find_if(match+1, seqs[k].end(), match_chr);
      ++i;
    }
  }

  // Overwrite original A3M alignment with new A2M alignment
  seqs = seqs_a2m;
}

template<class Alphabet>
void Alignment<Alphabet>::write(FILE* fout,
                                Format format,
                                int width) const {
  switch (format) {
    case FASTA:
    case A2M:
    case A3M:
      write_fasta_flavors(fout, format, width);
      break;
    case CLUSTAL:
    case PSI:
      write_clustal_flavors(fout, format, width);
      break;
    default:
      throw Exception("Unsupported alignment output format %i!", format);
  }
}

template<class Alphabet>
void Alignment<Alphabet>::write_fasta_flavors(FILE* fout,
                                              Format format,
                                              int width) const {
  for (int k = 0; k < num_seqs(); ++k) {
    fprintf(fout, ">%s\n", headers_[k].c_str());
    int j = 0;  // counts printed characters
    for (int i = 0; i < num_cols(); ++i) {
      switch (format) {
        case FASTA:
          fputc(to_match_chr(chr(k, i)), fout);
          ++j;
          break;
        case A2M:
          if (match_column_[i])
            fputc(to_match_chr(chr(k, i)), fout);
          else
            fputc(to_insert_chr(chr(k, i)), fout);
          ++j;
          break;
        case A3M:
          if (match_column_[i]) {
            fputc(to_match_chr(chr(k, i)), fout);
            ++j;
          } else if (seq(k, i) != Alphabet::instance().gap()
                     && seq(k, i) != Alphabet::instance().endgap()) {
            fputc(to_insert_chr(chr(k, i)), fout);
            ++j;
          }
          break;
        default:
          throw Exception("Unsupported alignment output format %i!", format);
      }
      if (j % width == 0) fputc('\n', fout);
    }
    if (j % width != 0) putc('\n', fout);
  }
}

template<class Alphabet>
void Alignment<Alphabet>::write_clustal_flavors(FILE* fout,
                                                Format format,
                                                int width) const {
  const int kHeaderWidth = 18;

  // Convert alignment to character representation
  std::vector<std::string> seqs(num_seqs(), "");
  for (int k = 0; k < num_seqs(); ++k) {
    for (int i = 0; i < num_cols(); ++i) {
      char c = Alphabet::instance().itoc(seqs_[i][k]);
      if (c != '-' && !match_column_[i] && format == PSI) c = to_insert_chr(c);
      seqs[k].push_back(c);
    }
  }

  // Print alignment in blocks
  if (format == CLUSTAL) fputs("CLUSTAL\n\n", fout);
  while (!seqs.front().empty()) {
    for (int k = 0; k < num_seqs(); ++k) {
      std::string header(headers_[k].substr(0, headers_[k].find_first_of(' ')));
      if (static_cast<int>(header.length()) <= kHeaderWidth) {
        header += std::string(kHeaderWidth - header.length() + 1, ' ');
        fputs(header.c_str(), fout);
      } else {
        fputs(header.substr(0, kHeaderWidth).c_str(), fout);
      }

      int len = std::min(width, static_cast<int>(seqs[k].length()));
      fputs(seqs[k].substr(0, len).c_str(), fout);
      fputc('\n', fout);
      seqs[k].erase(0, len);
    }
    fputc('\n', fout);  // blank line after each block
  }
}

template<class Alphabet>
void Alignment<Alphabet>::resize(int num_seqs, int num_cols) {
  if (num_seqs == 0 || num_cols == 0)
    throw Exception("Bad alignment dimensions: num_seqs=%i num_cols=%i",
                    num_seqs, num_cols);

  seqs_.resize(num_cols, num_seqs);
  column_indexes_.resize(num_cols);
  match_indexes.resize(num_cols);
  match_column_.resize(num_cols);
  headers_.resize(num_seqs);
}

template<class Alphabet>
void Alignment<Alphabet>::assign_match_columns_by_sequence(int k) {
  for (int i = 0; i < num_cols(); ++i)
    match_column_[i] = seqs_[i][k] < Alphabet::instance().gap();
  set_match_indexes();
}

template<class Alphabet>
void Alignment<Alphabet>::assign_match_columns_by_gap_rule(int gap_threshold) {
  LOG(DEBUG) << "Marking columns with more than " << gap_threshold
             << "% of gaps as insert columns ...";

  // Global weights are sufficient for calculation of gap percentage
  std::vector<float> wg;
  global_weights_and_diversity(*this, wg);
  for (int i = 0; i < num_cols(); ++i) {
    float gap = 0.0f;
    float res = 0.0f;

    for (int k = 0; k < num_seqs(); ++k)
      if (seqs_[i][k] < Alphabet::instance().any())
        res += wg[k];
      else
        gap += wg[k];

    float percent_gaps = 100.0f * gap / (res + gap);

    LOG(DEBUG3) << "percent gaps[" << i << "]=" << percent_gaps;

    match_column_[i] = percent_gaps <= static_cast<float>(gap_threshold);
  }
  set_match_indexes();

  LOG(DEBUG) << "num_cols=" << num_cols() << "  num_match_cols="
             << num_match_cols() << "  num_insert_cols=" << num_insert_cols();
}

template<class Alphabet>
void Alignment<Alphabet>::remove_insert_columns() {
  // Create new sequence matrix
  const int match_cols = num_match_cols();
  matrix<char> new_seqs(match_cols, num_seqs());
  for (int i = 0; i < match_cols; ++i) {
    for (int k = 0; k < num_seqs(); ++k) {
      new_seqs[i][k] = seqs_[match_indexes[i]][k];
    }
  }
  seqs_.resize(match_cols, num_seqs());
  seqs_ = new_seqs;

  // Update match indexes
  column_indexes_.resize(match_cols);
  match_indexes.resize(match_cols);
  match_column_.resize(match_cols);
  for (int i = 0; i < match_cols; ++i) {
    column_indexes_[i] = i;
    match_column_[i]   = true;
  }
  set_match_indexes();
}

template<class Alphabet>
void Alignment<Alphabet>::Merge(const Alignment<Alphabet>& ali) {
  if (num_match_cols() != ali.num_match_cols()) return;
  // FIXME: Keep insert columns when merging two alignments
  remove_insert_columns();

  // Copy and keep track of headers that are not already contained in master
  // alignment.
  std::vector<bool> include_seq(ali.num_seqs(), false);
  std::vector<std::string> headers_merged(headers_);
  for (int k = 0; k < ali.num_seqs(); ++k) {
    if (find(headers_.begin(), headers_.end(), ali.header(k)) == headers_.end()) {
      include_seq[k] = true;
      headers_merged.push_back(ali.header(k));
    }
  }

  // Append new sequences to alignment
  matrix<char> seqs_merged(num_match_cols(), headers_merged.size());
  for (int i = 0; i < num_match_cols(); ++i) {
    for (int k = 0; k < num_seqs(); ++k) {
      seqs_merged[i][k] = seqs_[match_indexes[i]][k];
    }
    int l = num_seqs(); // index of sequence k in merged alignment
    for (int k = 0; k < ali.num_seqs(); ++k) {
      if (include_seq[k]) {
        seqs_merged[i][l] = ali.seqs_[ali.match_indexes[i]][k];
        ++l;
      }
    }
  }

  // Apply the merging
  headers_ = headers_merged;
  seqs_    = seqs_merged;

  LOG(DEBUG) << "Merged alignment:";
  LOG(DEBUG) << *this;
}

template<class Alphabet>
Sequence<Alphabet> Alignment<Alphabet>::GetSequence(int k) const {
  std::string seq_str;
  for (int i = 0; i < num_cols(); ++i)
    if (seq(0,i) < Alphabet::instance().gap())
      seq_str += chr(0,i);

  return Sequence<Alphabet>(header(k), seq_str);
}

template<class Alphabet>
std::ostream& operator<< (std::ostream& out, const Alignment<Alphabet>& ali) {
  const int kHeaderWidth = 18;
  const int kLineLen     = 100;

  // Convert alignment to character representation
  std::vector<std::string> seqs(ali.num_seqs(), "");
  for (int k = 0; k < ali.num_seqs(); ++k) {
    for (int i = 0; i < ali.num_cols(); ++i) {
      char c = Alphabet::instance().itoc(ali.seqs_[i][k]);
      seqs[k].push_back(c);
    }
  }

  // Print alignment in blocks
  out << "CLUSTAL\n\n";
  while (!seqs.front().empty()) {
    for (int k = 0; k < ali.num_seqs(); ++k) {
      std::string header =
          ali.headers_[k].substr(0, ali.headers_[k].find_first_of(' '));
      if (static_cast<int>(header.length()) <= kHeaderWidth)
        out << header + std::string(kHeaderWidth - header.length(), ' ') << ' ';
      else
        out << header.substr(0, kHeaderWidth) << ' ';

      out << seqs[k].substr(0, std::min(kLineLen, static_cast<int>(seqs[k].length())))
          << std::endl;
      seqs[k].erase(0, std::min(kLineLen, static_cast<int>(seqs[k].length())));
    }
    out << std::endl;  // blank line after each block
  }
  return out;
}



template<class Alphabet>
float global_weights_and_diversity(const Alignment<Alphabet>& alignment,
                                   std::vector<float>& wg) {
  const float kZero = 1E-10;  // for calculation of entropy
  const int num_seqs = alignment.num_seqs();
  const int num_cols = alignment.num_match_cols();
  const int alphabet_size = Alphabet::instance().size();
  const int any   = Alphabet::instance().any();

  // Return values
  wg.assign(num_seqs, 0.0f);  // global sequence weights
  float neff = 0.0f;          // diversity of alignment

  std::vector<int> n(num_seqs, 0);          // Number of residues in sequence i
  std::vector<float> fj(alphabet_size, 0.0f);           // to calculate entropy
  std::vector<int> adiff(num_cols, 0);      // different letters in each column
  matrix<int> counts(num_cols, alphabet_size, 0);  // column counts (excl. ANY)

  LOG(INFO) << "Calculation of global weights and alignment diversity ...";

  // Count number of residues in each column
  for (int i = 0; i < num_cols; ++i)
    for (int k = 0; k < num_seqs; ++k)
      if (alignment[i][k] < any) {
        ++counts[i][alignment[i][k]];
        ++n[k];
      }

  // Count number of different residues in each column
  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      if (counts[i][a]) ++adiff[i];

  // Calculate weights
  for (int i = 0; i < num_cols; ++i)
    for (int k = 0; k < num_seqs; ++k)
      if (adiff[i] > 0 && alignment[i][k] < any)
        wg[k] += 1.0f/(adiff[i] * counts[i][alignment[i][k]] * n[k]);
  normalize_to_one(&wg[0], num_seqs);

  // Calculate number of effective sequences
  for (int i = 0; i < num_cols; ++i) {
    reset(&fj[0], alphabet_size);
    for (int k = 0; k < num_seqs; ++k)
      if (alignment[i][k] < any) fj[alignment[i][k]] += wg[k];
    normalize_to_one(&fj[0], alphabet_size);
    for (int a = 0; a < alphabet_size; ++a)
      if (fj[a] > kZero) neff -= fj[a] * fast_log2(fj[a]);
  }
  neff = pow(2.0, neff / num_cols);

  LOG(DEBUG) << "neff=" << neff;

  return neff;
}

template<class Alphabet>
std::vector<float> position_specific_weights_and_diversity(
    const Alignment<Alphabet>& alignment, matrix<float>& w) {
  // Maximal fraction of sequences with an endgap
  const float kMaxEndgapFraction = 0.1;
  // Minimum number of columns in subalignments
  const int kMinCols = 10;
  // Zero value for calculation of entropy
  const float kZero = 1E-10;
  const int num_seqs  = alignment.num_seqs();
  const int num_cols  = alignment.num_match_cols();
  const int alphabet_size  = Alphabet::instance().size();
  const int any    = Alphabet::instance().any();
  const int endgap = Alphabet::instance().endgap();

  // Return values
  std::vector<float> neff(num_cols, 0.0f);      // diversity of subalignment i
  w.resize(num_cols, num_seqs);                 // weight of seq k in column i
  w = matrix<float>(num_cols, num_seqs, 0.0f);  // init to zero

  // Helper variables
  int ncoli = 0;        // number of columns j that contribute to neff[i]
  int nseqi = 0;        // number of sequences in subalignment i
  int ndiff = 0;        // number of different alphabet letters
  bool change = false;  // has the set of sequences in subalignment changed?

  // Number of seqs with some residue in column i AND a at position j
  matrix<int> n(num_cols, endgap + 1, 0);
  // To calculate entropy
  std::vector<float> fj(alphabet_size, 0.0f);
  // Weight of sequence k in column i, calculated from subalignment i
  std::vector<float> wi(num_seqs, 0.0f);
  std::vector<float> wg;                       // global weights
  std::vector<int> nseqi_debug(num_cols, 0);   // debugging
  std::vector<int> ncoli_debug(num_cols, 0);   // debugging

  global_weights_and_diversity(alignment, wg);

  LOG(INFO) << "Calculation of position-specific weights ...";
  LOG(DEBUG2) << "i      ncol   nseq   neff";

  for (int i = 0; i < num_cols; ++i) {
    change = false;
    for (int k = 0; k < num_seqs; ++k) {
      if ((i == 0 || alignment[i-1][k] >= any) && alignment[i][k] < any) {
        change = true;
        ++nseqi;
        for (int j = 0; j < num_cols; ++j)
          ++n[j][alignment[j][k]];
      } else if (i > 0 && alignment[i-1][k] < any && alignment[i][k] >= any) {
        change = true;
        --nseqi;
        for (int j = 0; j < num_cols; ++j)
          --n[j][alignment[j][k]];
      }
    }  // for k over num_seqs
    nseqi_debug[i] = nseqi;

    if (change) {  // set of sequences in subalignment-inl.has changed
      ncoli = 0;
      reset(&wi[0], num_seqs);

      for (int j = 0; j < num_cols; ++j) {
        if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
        ndiff = 0;
        for (int a = 0; a < alphabet_size; ++a) if (n[j][a]) ++ndiff;
        if (ndiff == 0) continue;
        ++ncoli;
        for (int k = 0; k < num_seqs; ++k) {
          if (alignment[i][k] < any && alignment[j][k] < any) {
            if (n[j][alignment[j][k]] == 0) {
              LOG(WARNING) << "Mi=" << i << ": n[" << j << "][ali[" << j
                           << "][" << k << "]=0!";
            }
            wi[k] += 1.0f / static_cast<float>((n[j][alignment[j][k]] * ndiff));
          }
        }
      }  // for j over num_cols
      normalize_to_one(&wi[0], num_seqs);

      if (ncoli < kMinCols)  // number of columns in subalignment insufficient?
        for (int k = 0; k < num_seqs; ++k)
          if (alignment[i][k] < any)
            wi[k] = wg[k];
          else
            wi[k] = 0.0f;

      neff[i] = 0.0f;
      for (int j = 0; j < num_cols; ++j) {
        if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
        reset(&fj[0], alphabet_size);

        for (int k = 0; k < num_seqs; ++k)
          if (alignment[i][k] < any && alignment[j][k] < any)
            fj[alignment[j][k]] += wi[k];
        normalize_to_one(&fj[0], alphabet_size);
        for (int a = 0; a < alphabet_size; ++a)
          if (fj[a] > kZero) neff[i] -= fj[a] * fast_log2(fj[a]);
      }  // for j over num_cols

      neff[i] = (ncoli > 0 ? fast_pow2(neff[i] / ncoli) : 1.0f);

    } else {  // set of sequences in subalignment-inl.has NOT changed
      neff[i] = (i == 0 ? 0.0 : neff[i-1]);
    }

    for (int k = 0; k < num_seqs; ++k) w[i][k] = wi[k];
    ncoli_debug[i] = ncoli;

    LOG(DEBUG2) << strprintf("%7i%7i%7i%7.1f",
                             i, ncoli_debug[i], nseqi_debug[i], neff[i]);
  }  // for i over num_cols

  return neff;
}

template<class Alphabet>
typename Alignment<Alphabet>::Format alignment_format_from_string(
    const std::string& s) {
  if (s == "fas")
    return Alignment<Alphabet>::FASTA;
  else if (s == "a2m")
    return Alignment<Alphabet>::A2M;
  else if (s == "a3m")
    return Alignment<Alphabet>::A3M;
  else if (s == "clu")
    return Alignment<Alphabet>::CLUSTAL;
  else if (s == "psi")
    return Alignment<Alphabet>::PSI;
  else
    throw Exception("Unable to infer alignment format from extension '%s'!",
                    s.c_str());
}

}  // namespace cs

#endif  // SRC_ALIGNMENT_INL_H_
