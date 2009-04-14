// Copyright 2009, Andreas Biegert

#ifndef SRC_CONTEXT_SPECIFIC_BLAST_H_
#define SRC_CONTEXT_SPECIFIC_BLAST_H_

#include <map>
#include <memory>
#include <string>

#include "emitter.h"
#include "profile_library.h"
#include "library_pseudocounts.h"
#include "sequence-inl.h"

namespace cs {

// Options for CS-BLAST algorithm.
struct CSBlastOptions : public EmissionOptions {
  typedef std::map<char, std::string> psiblast_options;

  CSBlastOptions() { SetDefaults(); }
  virtual ~CSBlastOptions() { }

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

struct Pssm { }


// Implementation of CS-BLAST algorithm.
class CSBlast {
 public:
  // Constructor to compare a single sequence against a database of protein
  // sequences.
  CSBlast(const Sequence<AminoAcid>* query,
          const ProfileLibrary<AminoAcid>* lib,
          const CSBlastOptions* opts);

  // Constructor to compare a single sequence against a database of protein
  // sequences.
  CSBlast(const Sequence<AminoAcid>* query,
          const ProfileLibrary<AminoAcid>* lib,
          const CSBlastOptions* opts);

  ~ContextSpecificBlast() { }

  // Runs the CS-BLAST algorithm by adding context-specific pseudocounts to the
  // query sequence and then jumpstarting PSI-BLAST with the resulting profile.
  int Run();


  // Adds context-specific pseudocounts to alignment derived profile.
  virtual void add_to_profile(const Admixture& pca,
                              CountProfile<Alphabet>* p) const;

 private:
  // Profile library with context profiles.
  const ProfileLibrary<Alphabet>& lib_;
  // Needed to compute emission probabilities of context profiles.
  const Emitter<Alphabet> emitter_;

  DISALLOW_COPY_AND_ASSIGN(ContextSpecificBlast);
};  // ContextSpecificBlast

}  // namespace cs

#endif  // SRC_CONTEXT_SPECIFIC_BLAST_H_
