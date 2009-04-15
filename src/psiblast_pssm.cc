// Copyright 2009, Andreas Biegert

#include "psiblast_pssm.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "profile-inl.h"
#include "sequence-inl.h"

namespace cs {

PsiBlastPssm::PsiBlastPssm(Sequence<AminoAcid> query,
                           Profile<AminoAcid> profile)
    : query_(new Sequence<AminoAcid>(query)),
      profile_(new Profile<AminoAcid>(profile)) {}

PsiBlastPssm::PsiBlastPssm(FILE* fin) {
  Read(fin);
}

void PsiBlastPssm::Read(FILE* /* fin */) {
  // TODO
}

void PsiBlastPssm::Write(FILE* fout) {
  int query_length = query_->length();
  int bytes = 0;

  bytes = fwrite(reinterpret_cast<char*>(&query_length), 1, kIntSize, fout);
  for(int i = 0; i < query_length; ++i) {
    char c = AminoAcid::instance().itoc(query_->at(i));
    bytes = fwrite(&c, 1, kCharSize, fout);
  }

  for(int i = 0; i < query_length; ++i) {
    for(int a = 0; a < AminoAcid::instance().size(); ++a) {
      double p = (*profile_)[i][a];
      bytes = fwrite(reinterpret_cast<char*>(&p), 1, kCharSize, fout);
    }
  }
}

}  // namespace cs
