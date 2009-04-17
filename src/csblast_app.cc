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
#include "psiblast.h"
#include "psiblast_pssm.h"
#include "scoped_ptr.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;
using std::map;

namespace cs {

struct CSBlastAppOptions : public EmissionOptions {
  typedef map<char, string> PsiBlastOptions;

  CSBlastAppOptions() { SetDefaults(); }
  virtual ~CSBlastAppOptions() {}

  // Set csbuild default parameters
  void SetDefaults() {
    infile              = "";
    outfile             = "";
    libfile             = "";
    restartfile         = "";
    outformat           = 0;
    pc_admix            = 1.0f;
    pc_ali              = 10.0f;
    global_weights      = false;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (libfile.empty()) throw Exception("No profile library provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Library input file for restarting
  string libfile;
  // Checkpoint file for PSI-BLAST restarting
  string restartfile;
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
  // Query sequence
  scoped_ptr< Sequence<AminoAcid> > query_;
  // Profile library for pseudocounts
  scoped_ptr< ProfileLibrary<AminoAcid> > lib_;
  // Pseudocount factory
  scoped_ptr< LibraryPseudocounts<AminoAcid> > pc_;
  // PSI-BLAST engine
  scoped_ptr<PsiBlast> psiblast_;
  // PSSM for PSI-BLAST jumpstarting
  scoped_ptr<PsiBlastPssm> pssm_;
};  // class CSBlastApp



void CSBlastApp::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('R', "restart", opts_.restartfile, opts_.restartfile);
  *options >> Option('m', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('D', "context-pc", opts_.libfile, opts_.libfile);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

  // Put remaining arguments into PSI-BLAST options map
  for(GetOpt_pp::short_iterator it = options->begin(); it != options->end(); ++it) {
    putc(it.option(), stdout);
    putc('\n', stdout);
    if (!it.extracted())
      opts_.psiblast_opts[it.option()] = it.args().front();
  }


  opts_.Validate();
}

void CSBlastApp::print_description() const {
  fputs("Search with an amino acid sequence against protein databases for locally\n"
        "similar sequences.\n", stream());
}

void CSBlastApp::print_banner() const {
  fputs("Usage: csblast -i <infile> -D <library> [options]\n", stream());
}

void CSBlastApp::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Input file with query sequence");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Output file with search results (def=stdout)");
  fprintf(stream(), "  %-30s %s (def=off)\n", "-D, --context-pc <library>",
          "Path to profile library.");
  fprintf(stream(), "  %-30s %s\n", "-R, --restart <filename>",
          "Input file for CS-BLAST restart");
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-m, --outformat [0,11]",
          "Alignment view option", opts_.outformat);
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weighting");
}

int CSBlastApp::run() {
  // Read query sequence
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
  query_.reset(new Sequence<AminoAcid>(fin));
  fclose(fin);

  // Read profile library
  fin = fopen(opts_.libfile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.libfile.c_str());
  lib_.reset(new ProfileLibrary<AminoAcid>(fin));
  fclose(fin);

  // Setup context-specific pseudocounts generator
  pc_.reset(new LibraryPseudocounts<AminoAcid>(lib_.get(), opts_));

  // Read PSSM from restart file or init PSSM with query profile generated
  // with context-specific pseudocounts
  if (!opts_.restartfile.empty()) {
    fin = fopen(opts_.restartfile.c_str(), "rb");
    if (!fin) throw Exception("Unable to read file '%s'!",
                              opts_.restartfile.c_str());
    pssm_.reset(new PsiBlastPssm(fin));
    fclose(fin);
  } else {
    CountProfile<AminoAcid> profile(*query_);
    pc_->add_to_sequence(*query_, ConstantAdmixture(opts_.pc_admix), &profile);
    pssm_.reset(new PsiBlastPssm(query_->ToString(), profile));
  }

  // Setup PSI-BLAST engine
  psiblast_.reset(new PsiBlast(query_.get(), pssm_.get(),
                               opts_.psiblast_opts));

  // Run CS-BLAST until convergence or for maximum number of iterations
  psiblast_->Run();

  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSBlastApp().main(argc, argv, stdout, "cblast");
}
