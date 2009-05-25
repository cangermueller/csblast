// Copyright 2009, Andreas Biegert

#ifndef SRC_CSBLAST_H_
#define SRC_CSBLAST_H_

#include <cstdio>

#include <map>
#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "blast_hits.h"
#include "psiblast_pssm.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation of database searching with PSI-BLAST.
class CSBlast {
 public:
  typedef std::map<char, std::string> Options;

  // Constructor to compare a single sequence against a database of protein
  // sequences.
  CSBlast(const Sequence<AminoAcid>* query,
          const Options& opts);
  // Constructor for restarting CS-BLAST iterations with a previously generated
  // PSSM.
  CSBlast(const Sequence<AminoAcid>* query,
          const PsiBlastPssm* pssm,
          const Options& opts);
  ~CSBlast() {}

  // Runs one iteration of PSI-BLAST with cs-PSSM
  int Run(FILE* fout);
  // Runs one iteration of PSI-BLAST with cs-PSSM and returns found hits
  // in BlastHits object. Works for alignment format "-m 0" only!
  int Run(FILE* fout, BlastHits* hits = NULL);
  // Sets path to PSI-BLAST executable.
  void set_exec_path(std::string exec_path) {
    exec_path_ = exec_path;
    if (*exec_path.rbegin() != kDirSep) exec_path_ += kDirSep;
  }
  // Gets path to PSI-BLAST executable.
  std::string exec_path() const { return exec_path_; }
  // Sets position specific scoring matrix
  void set_pssm(const PsiBlastPssm* pssm) { pssm_ = pssm; }

 private:
  // Default path to PSI-BLAST executable
  static const char* kPsiBlastExec;
  // Reference string for output
  static const char* kCSBlastReference;
  // Options to be ignored for building command line string
  static const char* kIgnoreOptions;

  void WriteQuery(std::string filepath) const;
  void WriteCheckpoint(std::string filepath) const;
  std::string ComposeCommandString(std::string queryfile,
                                   std::string checkpointfile = "") const;

  // Query sequence provided by constructor
  const Sequence<AminoAcid>* query_;
  // PSI-BLAST position-specific scoring matrix
  const PsiBlastPssm* pssm_;
  // Options map with PSI-BLAST specific command-line arguments
  Options opts_;
  // Path to PSI-BLAST executable
  std::string exec_path_;

  DISALLOW_COPY_AND_ASSIGN(CSBlast);
};  // class CSBlast

}  // namespace cs

#endif  // SRC_CSBLAST_H_
