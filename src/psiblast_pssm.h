// Copyright 2009, Andreas Biegert

#ifndef SRC_PSIBLAST_PSSM_H_
#define SRC_PSIBLAST_PSSM_H_

#include <cstdio>
#include <cstdlib>

#include <string>

#include "amino_acid.h"
#include "profile.h"
#include "scoped_ptr.h"
#include "sequence.h"

namespace cs {

// Container class for PSI-BLAST PSSM and the associated query sequence.
class PsiBlastPssm {
 public:
  // Constructor to create a PSSM from the query and its sequence profile.
  PsiBlastPssm(Sequence<AminoAcid> query, Profile<AminoAcid> profile);
  // Constructor to create a PSSM from a PSI-BLAST checkpoint.
  PsiBlastPssm(FILE* fin);

  // Do nothing - scoped_ptr's do the job!
  ~PsiBlastPssm() {}

  // Gets query sequence
  const Sequence<AminoAcid>& query() const { return *query_; }
  // Gets sequence profile
  const Profile<AminoAcid>& profile() const { return *profile_; }
  // Overwrites profile with new profile.
  void set_profile(const Profile<AminoAcid> profile) {
    profile_.reset(new Profile<AminoAcid>(profile));
  }
  // Overwrites existing PSSM with PSSM from PSI-BLAST checkpoint.
  void Read(FILE* fin);
  // Writes PSSM in PSI-BLAST binary checkpoint format to stream.
  void Write(FILE* fout);

 private:
  // Query sequence with which search was started
  scoped_ptr< Sequence<AminoAcid> > query_;
  // Evolving sequence profile
  scoped_ptr< Profile<AminoAcid> > profile_;
};  // class PsiBlastPssm

}  // namespace cs

#endif  // SRC_PSIBLAST_PSSM_H_
