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
    libfile             = "";
    informat            = "auto";
    outformat           = "prf";
    pc_admix            = 1.0f;
    pc_ali              = 10.0f;
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
  // Library input file for restarting
  string libfile;
  // Input file format
  string informat;
  // Output file format
  string outformat;
  // Overall pseudocount admixture
  float pc_admix;
  // Constant in pseudocount calculation for alignments
  float pc_ali;
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
  virtual void parse_options(GetOpt_pp* options);
  // Prints options summary to stream.
  virtual void print_options() const;
  // Prints short application description.
  virtual void print_description() const;
  // Prints usage banner to stream.
  virtual void print_banner() const;
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
};  // class CSBuildApp



template<class Alphabet>
void CSBuildApp<Alphabet>::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('I', "informat", opts_.informat, opts_.informat);
  *options >> Option('O', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('M', "matchcol", opts_.matchcol_assignment,
                     opts_.matchcol_assignment);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('D', "context-pc", opts_.libfile, opts_.libfile);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

  opts_.Validate();

  if (opts_.outfile.empty())
    opts_.outfile = get_file_basename(opts_.infile, false) + "prf";
  if (opts_.informat == "auto")
    opts_.informat = get_file_ext(opts_.infile);
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_description() const {
  fputs("Build a profile, PSSM, or HMM from an alignment or sequence.\n",
        stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_banner() const {
  fputs("Usage: csbuild -i <infile> [options]\n", stream());
  fputs("       csbuild -i <infile> -D <library> [options]\n", stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_options() const {
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
  fprintf(stream(), "  %-30s %s (def=off)\n", "-D, --context-pc <library>",
          "Add context-specific pseudocounts with profile library");
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
  profile.write(fout);
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

  // Iniialize profile library if needed
  if (!opts_.libfile.empty()) {
    FILE* fin = fopen(opts_.libfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read from jumpstart file '%s'!",
                              opts_.libfile.c_str());
    lib_.reset(new ProfileLibrary<Alphabet>(fin));
    fclose(fin);
  }

  // Read input sequence/alignment
  if (opts_.informat == "seq") {  // build profile from sequence
    Sequence<Alphabet> seq(fin);
    CountProfile<Alphabet> profile(seq);

    if (lib_) {
      LibraryPseudocounts<Alphabet> pc(lib_.get(),
                                       opts_.weight_center,
                                       opts_.weight_decay);
      fprintf(stream(), "Adding context-specific pseudocounts (admix=%-.2f) ...\n",
              opts_.pc_admix);
      pc.add_to_sequence(seq, ConstantAdmixture(opts_.pc_admix), &profile);
    }
    if (opts_.outformat == "chk")
      WriteCheckpoint(seq, profile);
    else
      WriteProfile(profile);

  } else {  // build profile from alignment
    typename Alignment<Alphabet>::Format f =
      alignment_format_from_string<Alphabet>(opts_.informat);
    Alignment<Alphabet> ali(fin, f);
    if (f == Alignment<Alphabet>::FASTA) {
      if (opts_.matchcol_assignment == CSBuildAppOptions::kMatchColAssignByQuery)
        ali.assign_match_columns_by_sequence(0);
      else
        ali.assign_match_columns_by_gap_rule(opts_.matchcol_assignment);
    }
    CountProfile<Alphabet> profile(ali, !opts_.global_weights);

    if (lib_) {
      LibraryPseudocounts<Alphabet> pc(lib_.get(),
                                       opts_.weight_center,
                                       opts_.weight_decay);
      fprintf(stream(), "Adding context-specific pseudocounts (admix=%-.2f) ...\n",
              opts_.pc_admix);
      pc.add_to_profile(DivergenceDependentAdmixture(opts_.pc_admix,
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
