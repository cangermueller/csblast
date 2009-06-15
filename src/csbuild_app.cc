// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "alignment-inl.h"
#include "count_profile-inl.h"
#include "mult_emission.h"
#include "exception.h"
#include "getopt_pp.h"
#include "pseudocounts.h"
#include "hmm_pseudocounts-inl.h"
#include "library_pseudocounts-inl.h"
#include "profile_library-inl.h"
#include "psiblast_pssm.h"
#include "scoped_ptr.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSBuildAppOptions {
  static const int kMatchColAssignByQuery = -1;

  CSBuildAppOptions() { SetDefaults(); }
  virtual ~CSBuildAppOptions() {}

  // Set csbuild default parameters
  void SetDefaults() {
    infile              = "";
    outfile             = "";
    contextfile         = "";
    informat            = "auto";
    outformat           = "prf";
    pc_admix            = 1.0f;
    pc_ali              = 12.0f;
    pc_engine           = "auto";
    matchcol_assignment = kMatchColAssignByQuery;
    global_weights      = false;
    weight_center       = 1.3f;
    weight_decay        = 0.9f;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Input file with context profile library or HMM
  string contextfile;
  // Input file format
  string informat;
  // Output file format
  string outformat;
  // Overall pseudocount admixture
  float pc_admix;
  // Constant in pseudocount calculation for alignments
  float pc_ali;
  // Pseudocount engine
  string pc_engine;
  // Match column assignment for FASTA alignments
  int matchcol_assignment;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // Weight of central column in multinomial emission
  float weight_center;
  // Exponential decay of window weights
  float weight_decay;
};  // CSBuildAppOptions


template<class Alphabet>
class CSBuildApp : public Application {
 private:
  // Runs the csbuild application.
  virtual int Run();
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp* options);
  // Prints options summary to stream.
  virtual void PrintOptions() const;
  // Prints short application description.
  virtual void PrintBanner() const;
  // Prints usage banner to stream.
  virtual void PrintUsage() const;
  // Prints output format options.
  void PrintOutputFormatOptions() const;
  // Writes profile to outfile
  void WriteProfile(const CountProfile<Alphabet>& profile) const;
  // Writes PSI-BLAST checkpoint file
  void WriteCheckpoint(const Sequence<Alphabet>& query,
                       const CountProfile<Alphabet>& profile) const;

  // Parameter wrapper
  CSBuildAppOptions opts_;
  // Profile library for pseudocounts
  scoped_ptr< ProfileLibrary<Alphabet> > lib_;
  // HMM for pseudocounts
  scoped_ptr< HMM<Alphabet> > hmm_;
  // Pseudocount engine
  scoped_ptr< Pseudocounts<Alphabet> > pc_;
};  // class CSBuildApp



template<class Alphabet>
void CSBuildApp<Alphabet>::ParseOptions(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('I', "informat", opts_.informat, opts_.informat);
  *options >> Option('O', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('M', "matchcol", opts_.matchcol_assignment,
                     opts_.matchcol_assignment);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('D', "context-data", opts_.contextfile, opts_.contextfile);
  *options >> Option('p', "pc-engine", opts_.pc_engine, opts_.pc_engine);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

  opts_.Validate();

  if (opts_.outfile.empty())
    opts_.outfile = get_file_basename(opts_.infile, false) + "prf";
  if (opts_.informat == "auto")
    opts_.informat = get_file_ext(opts_.infile);
  if (opts_.pc_engine == "auto" && !opts_.contextfile.empty())
    opts_.pc_engine = get_file_ext(opts_.contextfile);
}

template<class Alphabet>
void CSBuildApp<Alphabet>::PrintBanner() const {
  fputs("Build a profile or PSSM from an alignment or sequence.\n",
        stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::PrintUsage() const {
  fputs("Usage: csbuild -i <infile> [options]\n", stream());
  fputs("       csbuild -i <infile> -D <context data> [options]\n", stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::PrintOptions() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <file>",
          "Input file with alignment or sequence");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with serialized profile (def: <basename>.prf)");
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-I, --informat seq|fas|...",
          "Input format: seq, fas, a2m, or a3m", opts_.informat.c_str());

  PrintOutputFormatOptions();

  fprintf(stream(), "  %-30s %s\n", "-M, --matchcol [0:100]",
          "Make all FASTA columns with less than X% gaps match columns");
  fprintf(stream(), "  %-30s %s\n", "", "(def: make columns with residue in "
          "first sequence match columns)");
  fprintf(stream(), "  %-30s %s (def=off)\n", "-D, --context-data <file>",
          "Add context-specific pseudocounts with profile library");
  // fprintf(stream(), "  %-30s %s (def=%s)\n", "-p, --pc-engine lib|hmm",
  //         "Specify engine for pseudocount generation", opts_.pc_engine.c_str());
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admixture for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weights (def=off)");
}

template<class Alphabet>
void CSBuildApp<Alphabet>::PrintOutputFormatOptions() const {
  /* There is only one output format for nucleotides! */
}

template<>
void CSBuildApp<AminoAcid>::PrintOutputFormatOptions() const {
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-O, --outformat prf|chk",
          "Output format: profile or PSI-BLAST checkpoint",
          opts_.outformat.c_str());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::WriteProfile(
    const CountProfile<Alphabet>& profile) const {
  FILE* fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Unable to write profile to output file '%s'!",
                    opts_.outfile.c_str());
  profile.Write(fout);
  fprintf(stream(), "Wrote profile with %i columns to %s\n",
          profile.num_cols(), opts_.outfile.c_str());
  fclose(fout);
}

template<class Alphabet>
void CSBuildApp<Alphabet>::WriteCheckpoint(
    const Sequence<Alphabet>& /* query */,
    const CountProfile<Alphabet>& /* profile */) const {
  /* do nothing */
}

template<>
void CSBuildApp<AminoAcid>::WriteCheckpoint(
    const Sequence<AminoAcid>& query,
    const CountProfile<AminoAcid>& profile) const {
  PsiBlastPssm pssm(query.ToString(), profile);

  FILE* fout = fopen(opts_.outfile.c_str(), "wb");
  if (!fout) throw Exception("Unable to write profile to output file '%s'!",
                             opts_.outfile.c_str());
  pssm.Write(fout);
  fprintf(stream(), "Wrote profile with %i columns as checkpoint file to %s\n",
          profile.num_cols(), opts_.outfile.c_str());
  fclose(fout);
}

template<class Alphabet>
int CSBuildApp<Alphabet>::Run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read from input file '%s'!",
                            opts_.infile.c_str());

  if (!opts_.contextfile.empty() && opts_.pc_engine == "lib") {
    // Iniialize profile library and library pseudocounts if needed
    FILE* fin = fopen(opts_.contextfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read from file '%s'!", opts_.contextfile.c_str());
    lib_.reset(new ProfileLibrary<Alphabet>(fin));
    fclose(fin);

    pc_.reset(new LibraryPseudocounts<Alphabet>(lib_.get(),
                                                opts_.weight_center,
                                                opts_.weight_decay));

  } else if (!opts_.contextfile.empty() && opts_.pc_engine == "hmm") {
    // Iniialize HMM and HMM pseudocounts if needed
    FILE* fin = fopen(opts_.contextfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read from file '%s'!", opts_.contextfile.c_str());
    hmm_.reset(new HMM<Alphabet>(fin));
    fclose(fin);

    pc_.reset(new HMMPseudocounts<Alphabet>(hmm_.get(),
                                            opts_.weight_center,
                                            opts_.weight_decay));
  }

  // Build profile from sequence
  if (opts_.informat == "seq") {
    Sequence<Alphabet> seq(fin);
    CountProfile<Alphabet> profile(seq);

    if (pc_) {
      fprintf(stream(), "Adding context-specific pseudocounts (admix=%-.2f) ...\n",
              opts_.pc_admix);
      pc_->AddPseudocountsToSequence(seq, ConstantAdmixture(opts_.pc_admix), &profile);
    }
    if (opts_.outformat == "chk")
      WriteCheckpoint(seq, profile);
    else
      WriteProfile(profile);

  } else {  // Build profile from alignment
    typename Alignment<Alphabet>::Format f =
      AlignmentFormatFromString<Alphabet>(opts_.informat);
    Alignment<Alphabet> ali(fin, f);
    if (f == Alignment<Alphabet>::FASTA) {
      if (opts_.matchcol_assignment == CSBuildAppOptions::kMatchColAssignByQuery)
        ali.AssignMatchColumnsBySequence(0);
      else
        ali.AssignMatchColumnsByGapRule(opts_.matchcol_assignment);
    }
    CountProfile<Alphabet> profile(ali, !opts_.global_weights);

    if (pc_) {
      fprintf(stream(), "Adding context-specific pseudocounts (admix=%-.2f) ...\n",
              opts_.pc_admix);
      pc_->AddPseudocountsToProfile(DivergenceDependentAdmixture(opts_.pc_admix,
                                                       opts_.pc_ali), &profile);
    }
    if (opts_.outformat == "chk")
      WriteCheckpoint(ali.GetSequence(0), profile);
    else
      WriteProfile(profile);
  }

  fclose(fin);
  return 0;
}

}  // namespace cs
