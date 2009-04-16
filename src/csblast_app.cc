// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <map>
#include <string>

#include "globals.h"
#include "amino_acid.h"
#include "alignment-inl.h"
#include "count_profile-inl.h"
#include "emitter.h"
#include "exception.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "profile_library-inl.h"
#include "psiblast_pssm.h"
#include "scoped_ptr.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;
using std::map;

namespace cs {

struct CSBlastAppOptions : public EmissionOptions {
  map<char, string> PsiBlastOptions;

  CSBlastAppOptions() { SetDefaults(); }
  virtual ~CSBlastAppOptions() {}

  // Set csbuild default parameters
  void SetDefaults() {
    infile              = "";
    outfile             = "";
    libfile             = "";
    outformat           = 0;
    pc_admix            = 1.0f;
    pc_ali              = 10.0f;
    global_weights      = false;
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
  // BLAST output format
  int outformat;
  // Overall pseudocount admixture
  float pc_admix;
  // Constant in pseudocount calculation for alignments
  float pc_ali;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // PSI-BLAST options map
  PsiBlastOptions psiblast_opts;
};  // struct CSBlastAppOptions


class CSBlastApp : public Application {
 private:


  // Runs the csbuild application.
  virtual int run();
  // Parses command line options.
  virtual void parse_options(GetOpt_pp* options);
  // Prints options summary to stream.
  virtual void print_options() const;
  // Prints short application description.
  virtual void print_description() const;
  // Prints usage banner to stream.
  virtual void print_banner() const;

  // Parameter wrapper
  CSBlastAppOptions opts_;
  // Profile library for pseudocounts
  scoped_ptr< ProfileLibrary<AminoAcid> > lib_;
  // PSI-BLAST engine
  scoped_ptr<PsiBlast> psiblast_;
};  // class CSBlastApp



void CSBuildApp::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('m', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('D', "context-pc", opts_.libfile, opts_.libfile);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

  opts_.Validate();
}

void CSBuildApp::print_description() const {
  fputs("Search with an amino acid sequence against protein databases for locally\n"
        "similar sequences.\n", stream());
}

void CSBuildApp::print_banner() const {
  fputs("Usage: csblast -i <infile> -D <library> [options]\n", stream());
}

void CSBuildApp::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Input file with alignment or sequence");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Output file with serialized profile (def: <basename>.prf)");
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-m, --outformat [0,11]",
          "Alignment view option identical to PSI-BLAST",
          opts_.outformat.c_str());
  fprintf(stream(), "  %-30s %s (def=off)\n", "-D, --context-pc <library>",
          "Add context-specific pseudocounts with profile library");
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admixture for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weighting");
}

int CSBuildApp::run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read from input file '%s'!",
                            opts_.infile.c_str());
  Sequence<AminoAcid> seq(fin);
  fclose(fin);

  FILE* fin = fopen(opts_.libfile.c_str(), "r");
  if (!fin) throw Exception("Unable to read from jumpstart file '%s'!",
                            opts_.libfile.c_str());
  lib_.reset(new ProfileLibrary<AminoAcid>(fin));
  fclose(fin);

  return 0;
}

}  // namespace cs
