// Copyright 2011, Andreas Biegert

#ifndef CS_CSLAST_H_
#define CS_CSLAST_H_

#include <map>

#include "last_pssm.h"
#include "sequence-inl.h"

namespace cs {

typedef std::map<char, std::string> CSLastOptions;

// Encapsulation of database searching with LAST.
class CSLast {
  public:
    // Constructor for running LAST with a context-specific PSSM
    CSLast(const std::string& dbfile,
           const LastPssm* pssm,
           const CSLastOptions& opts);

    // Runs one iteration of CS-LAST
    int Run(FILE* fout);

    // Sets path to PSI-LAST executable.
    void set_exec_path(std::string exec_path) {
        exec_path_ = exec_path;
        if (*exec_path.rbegin() != kDirSep) exec_path_ += kDirSep;
    }

    // Gets path to PSI-LAST executable.
    std::string exec_path() const { return exec_path_; }

    // Sets position specific scoring matrix
    void set_pssm(const LastPssm* pssm) { pssm_ = pssm; }

    // Sets command line options for PSI-LAST
    void set_options(const CSLastOptions& opts) { opts_ = opts; }

  private:
    // Default path to LAST executable
    static const char* kLastExec;
    // Options to be ignored for building command line string
    static const char* kIgnoreOptions;

    void WritePssm(std::string filepath) const;

    std::string ComposeCommandString(std::string pssmfile) const;

    // The database file
    std::string dbfile_;
    // LAST position-specific scoring matrix
    const LastPssm* pssm_;
    // Options map with LAST specific command-line arguments
    CSLastOptions opts_;
    // Path to LAST executable
    std::string exec_path_;

    DISALLOW_COPY_AND_ASSIGN(CSLast);
};  // class CSLast

}  // namespace cs

#endif  // CS_CSLAST_H_
