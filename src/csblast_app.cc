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
#include "csblast.h"
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
    outformat           = 0;
    pc_admix            = 0.95f;
    pc_ali              = 10.0f;
    global_weights      = false;
    weight_center       = 1.3;
    weight_decay        = 0.9;
    blast_path          = "";
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
  // Path to PSI-BLAST executable
  string blast_path;
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
  scoped_ptr<CSBlast> csblast_;
  // PSSM for PSI-BLAST jumpstarting
  scoped_ptr<PsiBlastPssm> pssm_;
};  // class CSBlastApp



void CSBlastApp::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  //  *options >> Option('m', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('D', "context-pc", opts_.libfile, opts_.libfile);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);
  *options >> Option(' ', "blast-path", opts_.blast_path, opts_.blast_path);
  *options >> Option(' ', "BLAST_PATH", opts_.blast_path, opts_.blast_path);

  // Put remaining arguments into PSI-BLAST options map
  for(GetOpt_pp::short_iterator it = options->begin(); it != options->end(); ++it) {
    if (!it.extracted())
      opts_.psiblast_opts[it.option()] = it.args().front();
  }

  opts_.Validate();
}

void CSBlastApp::print_description() const {
  fputs("Search with an amino acid sequence against protein databases for locally\n"
        "similar sequences.\n"
        "Biegert, A. and Soding, J. (2009), Sequence context-specific profiles for\n"
        "homology searching. Proc Natl Acad Sci USA, 106 (10), 3770-3775\n",
        stream());
}

void CSBlastApp::print_banner() const {
  fputs("Usage: csblast -i <infile> -D <profile library> [options] [blast options]\n",
        stream());
}

void CSBlastApp::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Input file with query sequence");
  fprintf(stream(), "  %-30s %s\n", "-D, --context-pc <library>",
          "Path to library with context profiles for cs pseudocounts");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Output file with search results (def=stdout)");
  fprintf(stream(), "  %-30s %s\n", "-d, --database <dbname>",
          "Protein database to search against (def=nr)");
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-m, --outformat [0,11]",
          "Alignment view option", opts_.outformat);
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weighting");
  fprintf(stream(), "  %-30s %s\n", "    --blast-path",
          "Set path to directory with PSI-BLAST executable");
}

int CSBlastApp::run() {
  int status = 0;

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

  // Construct PSSM with query profile generated with context-specific
  // pseudocounts
  CountProfile<AminoAcid> profile(*query_);
  pc_->add_to_sequence(*query_, ConstantAdmixture(opts_.pc_admix), &profile);
  pssm_.reset(new PsiBlastPssm(query_->ToString(), profile));

  // Setup PSI-BLAST engine
  csblast_.reset(new CSBlast(query_.get(),
                             pssm_.get(),
                             opts_.psiblast_opts));
  if (!opts_.blast_path.empty())
    csblast_->set_exec_path(opts_.blast_path);

  // Run one iteration of CS-BLAST
  if (opts_.outfile.empty()) {
    status = csblast_->Run(stream());
  } else {
    FILE* fout = fopen(opts_.outfile.c_str(), "w");
    status = csblast_->Run(fout);
    fclose(fout);
  }

  return status;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSBlastApp().main(argc, argv, stdout, "cblast");
}
