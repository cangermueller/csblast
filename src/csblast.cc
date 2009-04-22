// Copyright 2009, Andreas Biegert

#include "csblast.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "pseudocounts.h"
#include "sequence-inl.h"

using std::string;

namespace cs {

#ifdef _WIN32
const char* CSBlast::kPsiBlastExec = "blastpgp.exe";
#else
const char* CSBlast::kPsiBlastExec = "blastpgp";
#endif

const char* CSBlast::kIgnoreOptions = "ioR";

const char* CSBlast::kCSBlastReference =
  "Reference for sequence context-specific profiles:\n"
  "Biegert, Andreas and Soding, Johannes (2009), \n"
  "\"Sequence context-specific profiles for homology searching\", \n"
  "Proc Natl Acad Sci USA, 106 (10), 3770-3775.";

CSBlast::CSBlast(const Sequence<AminoAcid>* query,
                   const Options& opts)
    : query_(query), pssm_(NULL), opts_(opts), exec_path_() {}

CSBlast::CSBlast(const Sequence<AminoAcid>* query,
                   const PsiBlastPssm* pssm,
                   const Options& opts)
    : query_(query), pssm_(pssm), opts_(opts), exec_path_() {}

int CSBlast::Run(FILE* fout) {
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

  // Read PSI-BLAST output and print to stdout
  const int kLF = 0x0A;
  int c;
  bool print_reference =
    (opts_.find('m') == opts_.end() || opts_['m'] == "0") &&
    (opts_.find('T') == opts_.end() || opts_['T'] == "F");
  while (!feof(blast_out)) {
    if ((c = fgetc(blast_out)) != EOF && fout) {
      // Print CS-BLAST reference before BLAST reference
      if (print_reference && c == kLF) {
        fputs("\n\n", fout);
        fputs(kCSBlastReference, fout);
        print_reference = false;
      }

      fputc(c, fout);
      fflush(fout);
    }
  }
  int status = pclose(blast_out);

  // Cleanup temporary files
  remove(queryfile.c_str());
  if (pssm_) remove(checkpointfile.c_str());

  return status;
}

void CSBlast::WriteQuery(string filepath) const {
  FILE* fout = fopen(filepath.c_str(), "w");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  query_->write(fout);
  fclose(fout);
}

void CSBlast::WriteCheckpoint(string filepath) const {
  FILE* fout = fopen(filepath.c_str(), "wb");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  pssm_->Write(fout);
  fclose(fout);
}

string CSBlast::ComposeCommandString(string queryfile,
                                     string checkpointfile) const {
  string rv(exec_path_);
  rv = rv + kPsiBlastExec + " ";
  rv += " -i " + queryfile;
  if (pssm_) rv += " -R " + checkpointfile;

  for (Options::const_iterator it = opts_.begin(); it != opts_.end(); ++it) {
    if (!strchr(kIgnoreOptions, it->first))
      rv = rv + " -" + it->first + " " + it->second;
  }

  fputs(rv.c_str(), stderr);

  return rv;
}

}  // namespace cs
