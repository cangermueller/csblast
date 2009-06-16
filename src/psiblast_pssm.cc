// Copyright 2009, Andreas Biegert

#include "psiblast_pssm.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "profile-inl.h"
#include "sequence-inl.h"

using std::string;

namespace cs {

PsiBlastPssm::PsiBlastPssm (const std::string& query,
                            const Profile<AminoAcid>& profile)
    : query_(query),
      profile_(new Profile<AminoAcid>(profile)) {
  assert(profile.num_cols() == static_cast<int>(query.length()));
}

PsiBlastPssm::PsiBlastPssm(FILE* fin) {
  Read(fin);
}

void PsiBlastPssm::Read(FILE* fin) {
  const int alph_size = AminoAcid::instance().size();
  int count = 0;
  int query_length = 0;

  count += fread(reinterpret_cast<char*>(&query_length), kIntSize, 1, fin);
  assert(count == 1);
  assert(query_length > 0);

  query_.resize(query_length);
  profile_.reset(new Profile<AminoAcid>(query_length));

  for(int i = 0; i < query_length; ++i) {
    char c = '\0';
    count += fread(&c, kCharSize, 1, fin);
    assert(AminoAcid::instance().valid(c));
    query_[i] = c;
  }
  for(int i = 0; i < query_length; ++i) {
    for(int a = 0; a < alph_size; ++a) {
      double p = (*profile_)[i][a];
      count += fread(reinterpret_cast<char*>(&p), kDoubleSize, 1, fin);
      profile_->at(i,a) = p;
    }
  }
  assert(count == 1 + query_length + query_length * alph_size);
}

void PsiBlastPssm::Write(FILE* fout) const {
  LOG(INFO) << "Writing query and profile as PSI-BLAST checkpoint ...";
  LOG(INFO) << query_;
  LOG(INFO) << *profile_;

  const int alph_size = AminoAcid::instance().size();
  int query_length = query_.length();
  int count = 0;

  count += fwrite(reinterpret_cast<char*>(&query_length), kIntSize, 1, fout);
  for(int i = 0; i < query_length; ++i) {
    char c = query_[i];
    count += fwrite(&c, kCharSize, 1, fout);
  }

  for(int i = 0; i < query_length; ++i) {
    for(int a = 0; a < alph_size; ++a) {
      double p = profile_->at(i,a);
      count += fwrite(reinterpret_cast<char*>(&p), kDoubleSize, 1, fout);
    }
  }
  assert(count == 1 + query_length + query_length * alph_size);
}

}  // namespace cs
