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
#include "csblast.h"
#include "csblast_iteration.h"
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
  typedef map<char, string> PsiBlastOptions;

  CSBlastAppOptions() { SetDefaults(); }
  virtual ~CSBlastAppOptions() {}

  // Set csbuild default parameters
  void SetDefaults() {
    infile              = "";
    outfile             = "";
    libfile             = "";
    alifile             = "";
    checkpointfile      = "";
    outformat           = 0;
    pc_admix            = 0.95f;
    pc_ali              = 12.0f;
    global_weights      = false;
    weight_center       = 1.3;
    weight_decay        = 0.9;
    blast_path          = "";
    iterations          = 1;
    inclusion           = 0.002;
  }

  // Validates parameter settings and throws exception if needed.
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
  // Output file for multiple alignment of hits
  string alifile;
  // Output file for checkpointing
  string checkpointfile;
  // BLAST output format
  int outformat;
  // Overall pseudocount admixture
  float pc_admix;
  // Constant in pseudocount calculation for alignments
  float pc_ali;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // Path to PSI-BLAST executable
  string blast_path;
  // Maximum number of iterations to use in CSI-BLAST
  int iterations;
  // E-value threshold for inclusion in CSI-BLAST
  float inclusion;
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
  // Initializes all class members for CSI-BLAST searches
  void Init();
  // Writes current PSSM as PSI-BLAST checkpoint to checkpointfile
  void SavePssm() const;
  // Writes multiple alignment of hits to file
  void SaveAlignment() const;

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
  // Alignment of included sequences
  scoped_ptr< Alignment<AminoAcid> > ali_;
};  // class CSBlastApp



void CSBlastApp::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('C', "checkpoint", opts_.checkpointfile, opts_.checkpointfile);
  *options >> Option('m', "outformat", opts_.outformat, opts_.outformat);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  *options >> Option('D', "context-pc", opts_.libfile, opts_.libfile);
  *options >> Option('j', "iterations", opts_.iterations, opts_.iterations);
  *options >> Option('h', "inclusion", opts_.inclusion, opts_.inclusion);
  *options >> Option(' ', "alignhits", opts_.alifile, opts_.alifile);
  *options >> Option(' ', "blast-path", opts_.blast_path, opts_.blast_path);
  *options >> Option(' ', "BLAST_PATH", opts_.blast_path, opts_.blast_path);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

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
  fputs("Usage: csblast -i <infile> -D <library> [options] [blastpgp options]\n",
        stream());
}

void CSBlastApp::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <file>",
          "Input file with query sequence");
  fprintf(stream(), "  %-30s %s\n", "-D, --context-pc <file>",
          "Path to library with context profiles for cs pseudocounts");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with search results (def=stdout)");
  fprintf(stream(), "  %-30s %s\n", "-d, --database <dbname>",
          "Protein database to search against (def=nr)");
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-m, --outformat [0,11]",
          "Alignment view option", opts_.outformat);
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-j, --iterations [1,inf[",
          "Maximum number of iterations to use in  CSI-BLAST", opts_.iterations);
  fprintf(stream(), "  %-30s %s (def=%-.3f)\n", "-h, --inclusion [0,inf[",
          "E-value threshold for inclusion in  CSI-BLAST", opts_.inclusion);
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant for alignment pseudocounts in CSI-BLAST",
          opts_.pc_ali);
  fprintf(stream(), "  %-30s %s\n", "    --alignhits <file>",
          "Write FASTA multiple alignment of hits to file");
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weighting");
  fprintf(stream(), "  %-30s %s\n", "    --blast-path <path>",
          "Set path to directory with PSI-BLAST executable");
}

int CSBlastApp::run() {
  int status = 0;

  Init();
  CSBlastIteration itr(opts_.iterations);

  while (itr) {
    LOG(INFO) << strprintf("Starting iteration %i ...", itr.IterationNumber());

    SavePssm();

    // Run one iteration of CS-BLAST
    FILE* fout = opts_.outfile.empty() ? stream() : fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Unable to write file '%s'!", opts_.outfile.c_str());
    BlastHits hits;
    status = csblast_->Run(fout, &hits);
    if (!opts_.outfile.empty()) fclose(fout);
    if (status != 0 || hits.empty()) break;

    hits.Filter(opts_.inclusion);
    if (!hits.empty())
      ali_->Merge(Alignment<AminoAcid>(hits));
    LOG(INFO) << strprintf("Found %i sequences in iteration %i (E-value < %5.0E)",
                           hits.num_hits(), itr.IterationNumber(), opts_.inclusion);
    itr.Advance(hits);

    if (itr) {
      CountProfile<AminoAcid> ali_profile(*ali_, !opts_.global_weights);
      pc_->add_to_profile(DivergenceDependentAdmixture(opts_.pc_admix, opts_.pc_ali),
                         &ali_profile);
      pssm_.reset(new PsiBlastPssm(query_->ToString(), ali_profile));
      csblast_->set_pssm(pssm_.get());
    }
  }

  SaveAlignment();

  return status;
}

void CSBlastApp::Init() {
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

  // Setup PSSM of query profile with context-specific pseudocounts if no
  // restart file is provided
  if (opts_.psiblast_opts.find('R') == opts_.psiblast_opts.end()) {
    CountProfile<AminoAcid> profile(*query_);
    pc_->add_to_sequence(*query_, ConstantAdmixture(opts_.pc_admix), &profile);
    pssm_.reset(new PsiBlastPssm(query_->ToString(), profile));
  }

  // Setup CS-BLAST engine
  if (opts_.psiblast_opts.find('R') == opts_.psiblast_opts.end())
    csblast_.reset(new CSBlast(query_.get(), pssm_.get(), opts_.psiblast_opts));
  else
    csblast_.reset(new CSBlast(query_.get(), opts_.psiblast_opts));
  if (!opts_.blast_path.empty())
    csblast_->set_exec_path(opts_.blast_path);

  // Setup alignment of included sequences
  ali_.reset(new Alignment<AminoAcid>(*query_));
}

void CSBlastApp::SavePssm() const {
  if (!opts_.checkpointfile.empty() && pssm_) {
    FILE* fchk = fopen(opts_.checkpointfile.c_str(), "wb");
    if (!fchk)
      throw Exception("Unable to write file '%s'!", opts_.checkpointfile.c_str());
    pssm_->Write(fchk);
    fclose(fchk);
  }
}

void CSBlastApp::SaveAlignment() const {
  if (!opts_.alifile.empty()) {
    FILE* fali = fopen(opts_.alifile.c_str(), "w");
    if (!fali) throw Exception("Unable to write file '%s'!", opts_.alifile.c_str());
    ali_->write(fali, Alignment<AminoAcid>::FASTA);
    fclose(fali);
  }
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSBlastApp().main(argc, argv, stdout, "csblast");
}
