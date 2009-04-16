// Copyright 2009, Andreas Biegert

#include "csblast.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "pseudocounts.h"
#include "sequence-inl.h"

using std::string;

namespace cs {

#ifdef _WIN32
const char* PsiBlast::kPsiBlastExec = "blastpgp.exe";
#else
const char* PsiBlast::kPsiBlastExec = "blastpgp";
#endif

PsiBlast::PsiBlast(const Sequence<AminoAcid>* query,
                   const Options* opts)
    : query_(query), pssm_(NULL), opts_(opts), exec_path(kPsiBlastExec) {}

PsiBlast::PsiBlast(const Sequence<AminoAcid>* query,
                   const PsiBlastPssm& pssm,
                   const Options* opts)
    : query_(query), pssm_(pssm), opts_(opts), exec_path(kPsiBlastExec) {}

int PsiBlast::Run() {
  // First create unique basename for sequence and checkpoint file
  char name_template[] = "/tmp/csblast_XXXXXX";
  const int captured_fd = mkstemp(name_template);
  const string basename       = name_template;
  const string queryfile      = basename + ".seq";
  const string checkpointfile = basename + ".chk";

  // Create PSI-BLAST input files on the fly
  WriteQuery(queryfile);
  if (pssm_) WriteCheckpoint(checkpointfile);

  // Run PSI-BLAST with provided options
  FILE *blast_out;
  blast_out = popen("/cluster/bioprogs/blast32/blastpgp ?", "r");
  if (!blast_out) throw Exception("Error in execution of '%s'",
                                 "/cluster/bioprogs/blast32/blastpgp ?");
  char buffer[MB];
  while (!feof(blast_out)) {
    if (fgets(buffer, MB, blast_out))
      printf(buffer);
  }
  printf( "\nProcess returned %d\n", pclose(blast_out));
}

void PsiBlast::WriteQuery(string filepath) {
  FILE* fout = fopen(filepath.c_str(), "w");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  query_->write(fout);
  fclose(fout);
}

void PsiBlast::WriteCheckpoint(string filepath) {
  FILE* fout = fopen(filepath.c_str(), "wb");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  pssm_->Write(fout);
  fclose(fout);
}

}  // namespace cs
