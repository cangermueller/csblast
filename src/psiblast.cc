// Copyright 2009, Andreas Biegert

#include "psiblast.h"

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
//const char* PsiBlast::kPsiBlastExec = "blastpgp";
const char* PsiBlast::kPsiBlastExec = "/cluster/bioprogs/blast32/blastpgp";
#endif

PsiBlast::PsiBlast(const Sequence<AminoAcid>* query,
                   const Options& opts)
    : query_(query), pssm_(NULL), opts_(opts), exec_path_(kPsiBlastExec) {}

PsiBlast::PsiBlast(const Sequence<AminoAcid>* query,
                   const PsiBlastPssm* pssm,
                   const Options& opts)
    : query_(query), pssm_(pssm), opts_(opts), exec_path_(kPsiBlastExec) {}

int PsiBlast::Run() {
  // Create unique basename for sequence and checkpoint file
  char name_template[] = "/tmp/csblast_XXXXXX";
  const int captured_fd = mkstemp(name_template);
  if (!captured_fd) throw Exception("Unable to make unique basename!");
  const string basename       = name_template;
  const string queryfile      = basename + ".seq";
  const string checkpointfile = basename + ".chk";

  // Write PSI-BLAST input files
  WriteQuery(queryfile);
  if (pssm_) WriteCheckpoint(checkpointfile);

  // Run PSI-BLAST with provided options
  string command(ComposeCommandString(queryfile, checkpointfile));
  FILE *blast_out;
  blast_out = popen(command.c_str(), "r");
  if (!blast_out) throw Exception("Error executing '%s'", command.c_str());

  // Read PSI-BLAST output and print to stdout if no outfile given
  char buffer[MB];
  while (!feof(blast_out)) {
    if (fgets(buffer, MB, blast_out))
      printf("%s", buffer);
  }
  int status = pclose(blast_out);

  // Cleanup temporary files
  remove(queryfile.c_str());
  if (pssm_) remove(checkpointfile.c_str());

  return status;
}

void PsiBlast::WriteQuery(string filepath) const {
  FILE* fout = fopen(filepath.c_str(), "w");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  query_->write(fout);
  fclose(fout);
}

void PsiBlast::WriteCheckpoint(string filepath) const {
  FILE* fout = fopen(filepath.c_str(), "wb");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  pssm_->Write(fout);
  fclose(fout);
}

string  PsiBlast::ComposeCommandString(string queryfile,
                                       string checkpointfile) const {
  string rv(exec_path_);
  rv += " -i " + queryfile;
  if (pssm_) rv += " -R " + checkpointfile;

  for (Options::const_iterator it = opts_.begin(); it != opts_.end(); ++it)
    rv + rv + " -" + it->first + " " + it->second;

  puts(rv.c_str());

  return rv;
}

}  // namespace cs
