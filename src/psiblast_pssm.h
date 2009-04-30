// Copyright 2009, Andreas Biegert

#ifndef SRC_PSIBLAST_PSSM_H_
#define SRC_PSIBLAST_PSSM_H_

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <string>

#include "amino_acid.h"
#include "exception.h"
#include "profile.h"
#include "scoped_ptr.h"

namespace cs {

// Container class for PSI-BLAST PSSM and the associated query sequence.
class PsiBlastPssm {
 public:
  // Constructor to create a PSSM from query string and a sequence profile.
  PsiBlastPssm(const std::string& query,
               const Profile<AminoAcid>& profile);
  // Constructor to create a PSSM from a PSI-BLAST checkpoint file.
  PsiBlastPssm(FILE* fin);
  ~PsiBlastPssm() {}

  // Gets query sequence
  const std::string query() const { return query_; }
  // Gets sequence profile
  const Profile<AminoAcid>& profile() const { return *profile_; }
  // Overwrites query with new sequence string.
  void set_query(const std::string& seq) {
    assert(profile_->num_cols() == static_cast<int>(seq.length()));
    query_ = seq;
  }
  // Overwrites profile with new profile.
  void set_profile(const Profile<AminoAcid>& profile) {
    assert(profile.num_cols() == static_cast<int>(query_.length()));
    profile_.reset(new Profile<AminoAcid>(profile));
  }
  // Overwrites existing PSSM with PSSM from PSI-BLAST checkpoint.
  void Read(FILE* fin);
  // Writes PSSM in PSI-BLAST binary checkpoint format to stream.
  void Write(FILE* fout) const;

 private:
  // Query sequence with which search was started
  std::string query_;
  // Evolving sequence profile
  scoped_ptr< Profile<AminoAcid> > profile_;

  DISALLOW_COPY_AND_ASSIGN(PsiBlastPssm);
};  // class PsiBlastPssm

}  // namespace cs

#endif  // SRC_PSIBLAST_PSSM_H_
