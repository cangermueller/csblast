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

struct CSBlastAppOptions {
  typedef map<char, string> PsiBlastOptions;

  CSBlastAppOptions() { SetDefaults(); }
  virtual ~CSBlastAppOptions() {}

  void SetDefaults() {
    infile              = "";
    outfile             = "";
    libfile             = "";
    ali_infile          = "";
    ali_outfile         = "";
    checkpointfile      = "";
    pc_admix            = 0.95f;
    pc_ali              = 12.0f;
    global_weights      = false;
    weight_center       = 1.3;
    weight_decay        = 0.9;
    blast_path          = "";
    iterations          = 1;
    inclusion           = 0.002;
    best                = false;
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
  // Alignment file for starting from existing alignment
  string ali_infile;
  // Output file for multiple alignment of hits
  string ali_outfile;
  // Output file for checkpointing
  string checkpointfile;
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
  // Include only the best HSP per hit in alignment
  bool best;
  // Weight of central column in multinomial emission
  float weight_center;
  // Exponential decay of window weights
  float weight_decay;
  // PSI-BLAST options map
  PsiBlastOptions psiblast_opts;
};  // struct CSBlastAppOptions


class CSBlastApp : public Application {
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
  // Initializes all class members for CSI-BLAST searches
  void Init();
  // Writes current PSSM as PSI-BLAST checkpoint to checkpointfile
  void SavePssm() const;
  // Writes multiple alignment of hits to file
  void SaveAlignment() const;

  // Default number of one-line descriptions and alignments in BLAST output.
  // This should be large enough to ensure that the BLAST output parser can
  // extract all relevant sequences for inclusion in next CSI-BLAST iteration.
  static const int kNumOutputAlis = 5000;

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
  *options >> Option('B', "alifile", opts_.ali_infile, opts_.ali_infile);
  *options >> Option('C', "checkpoint", opts_.checkpointfile, opts_.checkpointfile);
  *options >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  *options >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  *options >> Option('D', "context-data", opts_.libfile, opts_.libfile);
  *options >> Option('j', "iterations", opts_.iterations, opts_.iterations);
  *options >> Option('h', "inclusion", opts_.inclusion, opts_.inclusion);
  *options >> Option(' ', "alignhits", opts_.ali_outfile, opts_.ali_outfile);
  *options >> Option(' ', "weight-center", opts_.weight_center,
                     opts_.weight_center);
  *options >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
  *options >> Option(' ', "blast-path", opts_.blast_path, opts_.blast_path);
  *options >> Option(' ', "BLAST_PATH", opts_.blast_path, opts_.blast_path);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);
  *options >> OptionPresent(' ', "best", opts_.best);

  // Put remaining arguments into PSI-BLAST options map
  for(GetOpt_pp::short_iterator it = options->begin(); it != options->end(); ++it) {
    if (!it.extracted())
      opts_.psiblast_opts[it.option()] = it.args().front();
  }

  // Set number of output one-line descriptions and alignments for CSI-BLAST
  if (opts_.iterations > 1) {
    if (opts_.psiblast_opts.find('v') == opts_.psiblast_opts.end() ||
        atoi(opts_.psiblast_opts['v'].c_str()) < kNumOutputAlis)
      opts_.psiblast_opts['v']  = strprintf("%i", kNumOutputAlis);
    if (opts_.psiblast_opts.find('b') == opts_.psiblast_opts.end() ||
        atoi(opts_.psiblast_opts['b'].c_str()) < kNumOutputAlis)
      opts_.psiblast_opts['b']  = strprintf("%i", kNumOutputAlis);
    opts_.psiblast_opts['m'] = "0";  // force -m 0 format needed for parsing
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
  fputs("Usage: csblast -i <infile> -D <context data> [options] [blastpgp options]\n",
        stream());
}

void CSBlastApp::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <file>",
          "Input file with query sequence");
  fprintf(stream(), "  %-30s %s\n", "-D, --context-data <file>",
          "Path to library with context profiles for cs pseudocounts");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with search results (def=stdout)");
  fprintf(stream(), "  %-30s %s\n", "-B, --alifile <file>",
          "Input alignment file for CSI-BLAST restart");
  fprintf(stream(), "  %-30s %s\n", "-d, --database <dbname>",
          "Protein database to search against (def=nr)");
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-j, --iterations [1,inf[",
          "Maximum number of iterations to use in CSI-BLAST", opts_.iterations);
  fprintf(stream(), "  %-30s %s (def=%-.3f)\n", "-h, --inclusion [0,inf[",
          "E-value threshold for inclusion in  CSI-BLAST", opts_.inclusion);
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(stream(), "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant for alignment pseudocounts in CSI-BLAST",
          opts_.pc_ali);
  fprintf(stream(), "  %-30s %s\n", "    --alignhits <file>",
          "Write multiple alignment of hits in PSI format to file");
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
          "Weight of central profile column", opts_.weight_center);
  fprintf(stream(), "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
          "Parameter for exponential decay of window weights", opts_.weight_decay);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weights (def=off)");
  fprintf(stream(), "  %-30s %s\n", "    --best",
          "Include only the best HSP per hit in alignment (def=off)");
  fprintf(stream(), "  %-30s %s\n", "    --blast-path <path>",
          "Set path to directory with PSI-BLAST executable");
}

int CSBlastApp::Run() {
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
    if (status != 0 || opts_.iterations == 1 || hits.empty()) break;

    hits.Filter(opts_.inclusion);
    if (!hits.empty() && !hits[0].hsps.empty())
      ali_->Merge(Alignment<AminoAcid>(hits, opts_.best));
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

  // Save alignments of hits if results were in -m0 format
  if (opts_.psiblast_opts.find('R') == opts_.psiblast_opts.end() ||
      opts_.psiblast_opts['m'] == "0")
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
  pc_.reset(new LibraryPseudocounts<AminoAcid>(lib_.get(),
                                               opts_.weight_center,
                                               opts_.weight_decay));

  // Setup PSSM of query profile with context-specific pseudocounts if no
  // restart file is provided
  if (opts_.psiblast_opts.find('R') == opts_.psiblast_opts.end()) {
    if (opts_.ali_infile.empty()) {
      CountProfile<AminoAcid> profile(*query_);
      pc_->add_to_sequence(*query_, ConstantAdmixture(opts_.pc_admix), &profile);
      pssm_.reset(new PsiBlastPssm(query_->ToString(), profile));

    } else {
      FILE* fin = fopen(opts_.ali_infile.c_str(), "r");
      Alignment<AminoAcid> input_alignment(fin, Alignment<AminoAcid>::PSI);
      input_alignment.AssignMatchColumnsBySequence(0);
      fclose(fin);

      CountProfile<AminoAcid> profile(input_alignment, !opts_.global_weights);
      pc_->add_to_profile(DivergenceDependentAdmixture(opts_.pc_admix, opts_.pc_ali),
                          &profile);
      pssm_.reset(new PsiBlastPssm(query_->ToString(), profile));
    }
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
  if (!opts_.ali_outfile.empty()) {
    FILE* fali = fopen(opts_.ali_outfile.c_str(), "w");
    if (!fali) throw Exception("Unable to write file '%s'!",
                               opts_.ali_outfile.c_str());
    ali_->Write(fali, Alignment<AminoAcid>::PSI);
    fclose(fali);
  }
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSBlastApp().main(argc, argv, stdout, "csblast");
}
