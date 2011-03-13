// Copyright 2009, Andreas Biegert

#ifndef CS_PSIBLAST_PSSM_H_
#define CS_PSIBLAST_PSSM_H_

#include "profile-inl.h"
#include "sequence-inl.h"

namespace cs {

// Container class for a PSI-BLAST position specific scoring matrix (PSSM) and its
// associated query sequence.
struct PsiBlastPssm {
  // Constructor to create a PSSM from query string and a sequence profile.
  PsiBlastPssm(const Sequence<AA>& seq, const Profile<AA>& prof)
      : query(seq), profile(prof) {}

  // Constructor to create a PSSM from a PSI-BLAST checkpoint file.
  PsiBlastPssm(FILE* fin);

  // Overwrites existing PSSM with PSSM from PSI-BLAST checkpoint.
  void Read(FILE* fin) {
    int count = 0, query_length = 0;

    count += fread(reinterpret_cast<char*>(&query_length), kIntSize, 1, fin);
    assert(count == 1);
    assert(query_length > 0);

    query.Resize(query_length);
    profile.Resize(query_length);

    for (int i = 0; i < query_length; ++i) {
      char c;
      count += fread(&c, kCharSize, 1, fin);
      assert(AA::kValidChar[c]);
      query[i] = c;
    }

    for(int i = 0; i < query_length; ++i) {
      for(size_t a = 0; a < AA::kSize; ++a) {
        double p;
        count += fread(reinterpret_cast<char*>(&p), kDoubleSize, 1, fin);
        profile[i][a] = p;
      }
    }
    assert(count == 1 + query_length + query_length * AA::kSize);
  }

  // Writes PSSM in PSI-BLAST binary checkpoint format to stream.
  void Write(FILE* fout) const {
    LOG(INFO) << "Writing PSI-BLAST checkpoint ...";
    LOG(INFO) << query;
    LOG(INFO) << profile;

    int query_length = query.length();
    int count = 0;

    count += fwrite(reinterpret_cast<char*>(&query_length), kIntSize, 1, fout);
    for(int i = 0; i < query_length; ++i) {
      char c = query_[i];
      count += fwrite(&c, kCharSize, 1, fout);
    }

    for(int i = 0; i < query_length; ++i) {
      for(size_t a = 0; a < AA::kSize; ++a) {
        double p = profile[i][a];
        count += fwrite(reinterpret_cast<char*>(&p), kDoubleSize, 1, fout);
      }
    }
    assert(count == 1 + query_length + query_length * alph_size);
  }

  // Query sequence with which search was started
  Sequence<AA> query;
  // Evolving sequence profile
  Profile<AA> profile;
};  // class PsiBlastPssm

}  // namespace cs

#endif  // CS_PSIBLAST_PSSM_H_
