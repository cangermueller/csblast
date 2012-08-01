/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_CSBLAST_H_
#define CS_CSBLAST_H_

#include <map>

#include "blast_hits.h"
#include "pssm.h"
#include "sequence-inl.h"

namespace cs {

typedef std::map<char, std::string> CSBlastOptions;

// Encapsulation of database searching with PSI-BLAST.
class CSBlast {
 public:
  // Constructor to compare a single sequence against a database of protein
  // sequences.
  CSBlast(const Sequence<AA>* query,
          const CSBlastOptions& opts);

  // Constructor for restarting CS-BLAST iterations with a previously generated
  // PSSM.
  CSBlast(const Sequence<AA>* query,
          const Pssm* pssm,
          const CSBlastOptions& opts);

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
  void set_pssm(const Pssm* pssm) { pssm_ = pssm; }

  // Sets command line options for PSI-BLAST
  void set_options(const CSBlastOptions& opts) { opts_ = opts; }

  // Sets BLAST call emulation
  void set_emulate(bool emulate = true) { emulate_ = emulate; }

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
  const Sequence<AA>* query_;
  // PSI-BLAST position-specific scoring matrix
  const Pssm* pssm_;
  // Options map with PSI-BLAST specific command-line arguments
  CSBlastOptions opts_;
  // Path to PSI-BLAST executable
  std::string exec_path_;
  // Emulate BLAST call
  bool emulate_;

  DISALLOW_COPY_AND_ASSIGN(CSBlast);
};  // class CSBlast

}  // namespace cs

#endif  // CS_CSBLAST_H_
