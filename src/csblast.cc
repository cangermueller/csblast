// Copyright 2009, Andreas Biegert

#include "csblast.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "blast_hits.h"
#include "log.h"
#include "pseudocounts.h"
#include "sequence-inl.h"

using std::string;

namespace cs {

#ifdef _WIN32
const char* CSBlast::kPsiBlastExec = "blastpgp.exe";
#else
const char* CSBlast::kPsiBlastExec = "blastpgp";
#endif

const char* CSBlast::kIgnoreOptions = "o";

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

int CSBlast::Run(FILE* fout, BlastHits* hits) {
  // Create unique basename for sequence and checkpoint file
  char name_template[] = "/tmp/csblast_XXXXXX";
  const int captured_fd = mkstemp(name_template);
  if (!captured_fd) throw Exception("Unable to create unique filename!");
  const string basename       = name_template;
  const string queryfile      = basename + ".seq";
  const string checkpointfile = basename + ".chk";
  const string resultsfile    = basename + ".out";

  WriteQuery(queryfile);
  WriteCheckpoint(checkpointfile);
  FILE* fres = fopen(resultsfile.c_str(),"w+");
  if (!fres) throw Exception("Unable to open file '%s'!", resultsfile.c_str());

  // Run PSI-BLAST with provided options
  string command(ComposeCommandString(queryfile, checkpointfile));
  FILE* blast_out = popen(command.c_str(), "r");
  if (!blast_out) throw Exception("Error executing '%s'", command.c_str());

  // Read PSI-BLAST output and print to stdout
  const int kLF = 0x0A;
  int c;
  bool print_reference = ((opts_.find('m') == opts_.end() || opts_['m'] == "0") &&
                          (opts_.find('T') == opts_.end() || opts_['T'] == "F"));
  while (!feof(blast_out)) {
    if ((c = fgetc(blast_out)) != EOF && fout) {
      if (print_reference && c == kLF) {
        fputs("\n\n", fout);
        fputs(kCSBlastReference, fout);
        print_reference = false;
      }

      fputc(c, fout);
      fflush(fout);
      fputc(c, fres);
      fflush(fres);
    }
  }
  int status = pclose(blast_out);

  // Parse hits from PSI-BLAST results
  rewind(fres);
  hits->Read(fres);
  fclose(fres);

  // Cleanup temporary files
  remove(queryfile.c_str());
  remove(resultsfile.c_str());
  remove(checkpointfile.c_str());

  return status;
}

void CSBlast::WriteQuery(string filepath) const {
  FILE* fout = fopen(filepath.c_str(), "w");
  if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
  query_->Write(fout);
  fclose(fout);
}

void CSBlast::WriteCheckpoint(string filepath) const {
  if (pssm_) {
    FILE* fout = fopen(filepath.c_str(), "wb");
    if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
    pssm_->Write(fout);
    fclose(fout);
  }
}

string CSBlast::ComposeCommandString(string queryfile, string checkpointfile) const {
  string rv(exec_path_ + kPsiBlastExec);
  Options full_opts(opts_);
  full_opts['i'] = queryfile;
  if (pssm_) full_opts['R'] = checkpointfile;

  for (Options::const_iterator it = full_opts.begin(); it != full_opts.end(); ++it) {
    if (!strchr(kIgnoreOptions, it->first))
      rv = rv + " -" + it->first + " " + it->second;
  }
  LOG(INFO) << rv;

  return rv;
}

}  // namespace cs
