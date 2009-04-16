// Copyright 2009, Andreas Biegert

#ifndef SRC_CSBLAST_H_
#define SRC_CSBLAST_H_

#include <map>
#include <string>

#include "amino_acid.h"
#include "emitter.h"
#include "profile_library.h"
#include "psiblast_pssm.h"
#include "library_pseudocounts.h"
#include "sequence-inl.h"

namespace cs {

// Options for CS-BLAST algorithm.
struct CSBlastOptions : public EmissionOptions {
  typedef std::map<char, std::string> psiblast_options;

  CSBlastOptions() { SetDefaults(); }
  virtual ~CSBlastOptions() {}

  // Set CS-BLAST default parameters
  void SetDefaults() {
    pc_admix       = 1.0f;
    pc_ali         = 10.0f;
    global_weights = false;
    working_dir    = "";
#ifdef _WIN32
    psiblast_exec  = "blastpgp.exe"
#else
    psiblast_exec  = "blastpgp"
#endif
  }

  // Overall pseudocount admixture
  float pc_admix;
  // Constant in pseudocount calculation for alignments
  float pc_ali;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // Working directory for temporary files
  std::string working_dir;
  // PSI-BLAST executable or absolute path to it
  std::string psiblast_exec;
};

// Implementation of CS-BLAST algorithm.
class CSBlast {
 public:
  // Constructor to compare a single sequence against a database of protein
  // sequences.
  CSBlast(const Sequence<AminoAcid>& query,
          const Pseudocounts<AminoAcid>* pc,
          const CSBlastOptions* opts);
  // Constructor for restarting CS-BLAST iterations with a previously generated
  // PSSM.
  CSBlast(const PsiBlastPssm& pssm,
          const Pseudocounts<AminoAcid>* pc,
          const CSBlastOptions* opts);
  ~CSBlast() {}

  // Runs the CS-BLAST algorithm by adding context-specific pseudocounts to the
  // query sequence and then jumpstarting PSI-BLAST with the resulting profile.
  int Run();

 private:
  // Query sequence (either extracted from PSSM or provided in constructor)
  const Sequence<AminoAcid> query_;
  // PSI-BLAST position-specific scoring matrix
  PsiBlatPssm pssm_;
  // Pseudocount factory for generation of sequence profile.
  const Pseudocounts<AminoAcid>* pc_;

  DISALLOW_COPY_AND_ASSIGN(CSBlast);
};  // ContextSpecificBlast

}  // namespace cs

#endif  // SRC_CSBLAST_H_
