// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "csblast.h"
#include "blast_hits.h"

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

CSBlast::CSBlast(const Sequence<AA>* query,
                 const CSBlastOptions& opts)
    : query_(query), pssm_(NULL), opts_(opts), exec_path_() {}

CSBlast::CSBlast(const Sequence<AA>* query,
                 const Pssm* pssm,
                 const CSBlastOptions& opts)
    : query_(query), pssm_(pssm), opts_(opts), exec_path_() {}

int CSBlast::Run(FILE* fout, BlastHits* hits) {
  int status = 0;
  string basename, queryfile, checkpointfile, resultsfile;

  try {
    // Create unique basename for sequence and checkpoint file
    char name_template[] = "/tmp/csblast_XXXXXX";
    const int captured_fd = mkstemp(name_template);
    if (!captured_fd) throw Exception("Unable to create unique filename!");
    basename       = name_template;
    queryfile      = basename + ".seq";
    checkpointfile = basename + ".chk";
    resultsfile    = basename + ".out";

    WriteQuery(queryfile);
    WriteCheckpoint(checkpointfile);

    // Run PSI-BLAST with provided options
    string command(ComposeCommandString(queryfile, checkpointfile));
    if (emulate_) {
      fprintf(fout, "%s\n", command.c_str());
      if (!basename.empty()) remove(basename.c_str());
      if (!queryfile.empty()) remove(queryfile.c_str());
      if (!checkpointfile.empty()) remove(checkpointfile.c_str());
      return 0;
    }

    FILE* fres = fopen(resultsfile.c_str(),"w+");
    if (!fres) throw Exception("Unable to open file '%s'!", resultsfile.c_str());
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
    status = pclose(blast_out);

    // Parse hits from PSI-BLAST results
    if (hits && (opts_.find('m') == opts_.end() || opts_['m'] == "0")) {
      rewind(fres);
      hits->Read(fres);
      fclose(fres);
    }

    // Cleanup temporary files
    if (!basename.empty()) remove(basename.c_str());
    if (!queryfile.empty()) remove(queryfile.c_str());
    if (!checkpointfile.empty()) remove(checkpointfile.c_str());
    if (!resultsfile.empty()) remove(resultsfile.c_str());

  } catch(const std::exception&) { // something went wrong => cleanup before rethrow
    // Cleanup temporary files
    if (!basename.empty()) remove(basename.c_str());
    if (!queryfile.empty()) remove(queryfile.c_str());
    if (!checkpointfile.empty()) remove(checkpointfile.c_str());
    if (!resultsfile.empty()) remove(resultsfile.c_str());
    throw;
  }
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

string CSBlast::ComposeCommandString(string queryfile,
                                     string checkpointfile) const {
  string rv(exec_path_ + kPsiBlastExec);
  if (!IsRegularFile(rv))
    throw Exception("No PSI-BLAST binary in directory '%s'!", exec_path_.c_str());
  CSBlastOptions full_opts(opts_);
  full_opts['i'] = queryfile;
  if (pssm_) full_opts['R'] = checkpointfile;

  typedef CSBlastOptions::const_iterator OptsIter;
  for (OptsIter it = full_opts.begin(); it != full_opts.end(); ++it)
    if (it->first == 'd')
      rv = rv + " -d '" + it->second + "'";
    else if (!strchr(kIgnoreOptions, it->first))
      rv = rv + " -" + it->first + " " + it->second;
  LOG(INFO) << rv;

  return rv;
}

}  // namespace cs
