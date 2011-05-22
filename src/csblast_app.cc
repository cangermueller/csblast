// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "alignment-inl.h"
#include "application.h"
#include "blast_hits.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf_pseudocounts-inl.h"
#include "crf-inl.h"
#include "lr_pseudocounts-inl.h"
#include "count_profile-inl.h"
#include "csblast.h"
#include "csblast_iteration.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "pssm.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;
using std::map;
using std::vector;

namespace cs {

typedef vector<Sequence<AA> > SeqVec;

struct CSBlastAppOptions {

  CSBlastAppOptions() { Init(); }

  void Init() {
    pc_admix        = 0.90;
    pc_ali          = 12.0;
    pc_neff         = 0.0;
    pc_engine       = "auto";
    global_weights  = false;
    weight_center   = 1.6;
    weight_decay    = 0.85;
    iterations      = 1;
    inclusion       = 0.002;
    best            = false;
    ndescr          = 500;
    nalis           = 250;
    emulate         = false;
    // penalty_alpha       = 0.0;
    // penalty_beta        = 0.1;
    // penalty_score_min   = 8.0;
    // no_penalty          = false;
  }

  // Validates parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (modelfile.empty()) throw Exception("No context data provided!");
    if (pc_admix <= 0 || pc_admix > 1.0) throw Exception("Pseudocounts admix invalid!");
    if (pc_neff < 1.0 && pc_neff != 0.0) 
      throw Exception("Target Neff for pseudocounts admixture invalid!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Input file with context profile library or HMM
  string modelfile;
  // Alignment file for starting from existing alignment
  string ali_infile;
  // Output file for multiple alignment of hits
  string ali_outfile;
  // Output file for checkpointing
  string checkpointfile;
  // Overall pseudocount admixture
  double pc_admix;
  // Constant in pseudocount calculation for alignments
  double pc_ali;
  // Target Neff for pseudocounts admixture
  double pc_neff;
  // Pseudocount engine
  string pc_engine;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // Path to PSI-BLAST executable
  string blast_path;
  // Maximum number of iterations to use in CSI-BLAST
  int iterations;
  // E-value threshold for inclusion in CSI-BLAST
  double inclusion;
  // Include only the best HSP per hit in alignment
  bool best;
  // Weight of central column in multinomial emission
  double weight_center;
  // Exponential decay of window weights
  double weight_decay;
  // Number of database sequences to show one-line descriptions for.
  int ndescr;
  // Number of database sequence to show alignments.
  int nalis;
  // Emulate BLAST call
  bool emulate;
  // Baseline penalty for adjusting E-values.
  // double penalty_alpha;
  // Repeat penalty strength
  // double penalty_beta;
  // Minimal score for repeat penalty to take effect
  // double penalty_score_min;
  // Perform score penalty for repeats
  // bool no_penalty;
  // PSI-BLAST options map
  CSBlastOptions csblast;
};  // struct CSBlastAppOptions


class CSBlastApp : public Application {
 private:
  // Runs the csbuild application.
  virtual int Run();
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& ops);
  // Prints options summary to stream.
  virtual void PrintOptions() const;
  // Prints short application description.
  virtual void PrintBanner() const;
  // Prints usage banner to stream.
  virtual void PrintUsage() const;
  // Initializes all class members for CSI-BLAST searches
  void Init();
  // Writes current PSSM as PSI-BLAST checkpoint to checkpointfile
  void SavePssm() const;
  // Writes multiple alignment of hits to file
  void SaveAlignment() const;
  // Add pseudocounts and prepares CS-BLAST engine for run with given query
  void PrepareForRun(const Sequence<AA>& query);

  // Default number of one-line descriptions and alignments in BLAST output.
  // This should be large enough to ensure that the BLAST output parser can
  // extract all relevant sequences for inclusion in next CSI-BLAST iteration.
  static const int kNumOutputAlis = 5000;

  // Parameter wrapper
  CSBlastAppOptions opts_;
  // Query sequence
  scoped_ptr<Sequence<AA> > query_;
  // Profile library for pseudocounts
  scoped_ptr<ContextLibrary<AA> > lib_;
  // CRF for pseudocounts
  scoped_ptr<Crf<AA> > crf_;
  // lr params for pseudocounts
  scoped_ptr<LrParams<AA> > lrparams_;
  // Pseudocount engine
  scoped_ptr<Pseudocounts<AA> > pc_;
  // PSI-BLAST engine
  scoped_ptr<CSBlast> csblast_;
  // PSSM for PSI-BLAST jumpstarting
  scoped_ptr<Pssm> pssm_;
  // Alignment of included sequences
  scoped_ptr<Alignment<AA> > ali_;
  // Vector with pointers to query sequences
  SeqVec queries_;
  // Repeat penalizer
  // scoped_ptr<RepeatPenalizer<AA> > penalizer_;
};  // class CSBlastApp



void CSBlastApp::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('B', "alifile", opts_.ali_infile, opts_.ali_infile);
  ops >> Option('C', "checkpoint", opts_.checkpointfile, opts_.checkpointfile);
  ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  ops >> Option('z', "pc-neff", opts_.pc_neff, opts_.pc_neff);
  ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
  ops >> Option('j', "iters", opts_.iterations, opts_.iterations);
  ops >> Option('h', "incl-evalue", opts_.inclusion, opts_.inclusion);
  ops >> Option('v', "descr", opts_.ndescr, opts_.ndescr);
  ops >> Option('b', "alis", opts_.nalis, opts_.nalis);
  ops >> Option(' ', "pc-engine", opts_.pc_engine, opts_.pc_engine);
  ops >> Option(' ', "alignhits", opts_.ali_outfile, opts_.ali_outfile);
  ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
  ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
  ops >> Option(' ', "blast-path", opts_.blast_path, opts_.blast_path);
  ops >> Option(' ', "BLAST_PATH", opts_.blast_path, opts_.blast_path);
  ops >> OptionPresent(' ', "emulate", opts_.emulate);
  ops >> OptionPresent(' ', "global-weights", opts_.global_weights);
  ops >> OptionPresent(' ', "best", opts_.best);

  // Put remaining arguments into PSI-BLAST options map
  for(GetOpt_pp::short_iterator it = ops.begin(); it != ops.end(); ++it) {
    if (!it.extracted())
      opts_.csblast[it.option()] = it.args().front();
  }

  if (opts_.iterations > 1)
    opts_.csblast['m'] = "0";  // force -m 0 format needed for parsing

  opts_.Validate();

  if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
    opts_.pc_engine = GetFileExt(opts_.modelfile);
}

void CSBlastApp::PrintBanner() const {
  fputs("Search with an amino acid sequence against protein databases for locally\n"
        "similar sequences.\n"
        "Biegert, A. and Soding, J. (2009), Sequence context-specific profiles\n"
        "for homology searching. Proc Natl Acad Sci USA, 106 (10), 3770-3775\n",
        out_);
}

void CSBlastApp::PrintUsage() const {
  fputs("Usage: csblast -i <queryfile> -D <contextdata> --blast-path <psiblastdir>"
        " [options] [blastpgp options]\n", out_);
}

void CSBlastApp::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with query sequence");
  fprintf(out_, "  %-30s %s\n", "-D, --context-data <file>",
          "Path to profile library with context profiles");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with search results (def=stdout)");
  fprintf(out_, "  %-30s %s\n", "-B, --alifile <file>",
          "Input alignment file for CSI-BLAST restart");
  fprintf(out_, "  %-30s %s\n", "-d, --database <dbname>",
          "Protein database to search against (def=nr)");
  fprintf(out_, "  %-30s %s (def=%i)\n", "-j, --iters [1,inf[",
          "Maximum number of iterations to use in CSI-BLAST", opts_.iterations);
  fprintf(out_, "  %-30s %s (def=%-.4f)\n", "-h, --incl-evalue [0,inf[",
          "E-value threshold for inclusion in  CSI-BLAST", opts_.inclusion);
  fprintf(out_, "  %-30s %s (def=%i)\n", "-v, --descr [1,inf[",
          "Number of sequences to show one-line descriptions for", opts_.ndescr);
  fprintf(out_, "  %-30s %s (def=%i)\n", "-b, --alis [1,inf[",
          "Number of sequences to show alignments for", opts_.nalis);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix ]0,1]",
          "Pseudocount admix for context-specific pseudocounts", opts_.pc_admix);
  fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant for alignment pseudocounts in CSI-BLAST", opts_.pc_ali);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-z, --pc-neff [1,inf[",
          "Target Neff for pseudocounts admixture", opts_.pc_neff);
  fprintf(out_, "  %-30s %s (def=%s)\n", "    --pc-engine <engine>",
           "Specify engine for pseudocount generation", opts_.pc_engine.c_str());
  fprintf(out_, "  %-30s %s\n", "", "<engine> = auto|crf|lib|lr");
  fprintf(out_, "  %-30s %s\n", "    --alignhits <file>",
          "Write multiple alignment of hits in PSI format to file");
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
          "Weight of central profile column", opts_.weight_center);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
          "Parameter for exponential decay of window weights", opts_.weight_decay);
  fprintf(out_, "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weights (def=off)");
  fprintf(out_, "  %-30s %s\n", "    --best",
          "Include only the best HSP per hit in alignment (def=off)");
  fprintf(out_, "  %-30s %s\n", "    --emulate",
          "Emulate BLAST call (def=off)");
  // fprintf(out_, "  %-30s %s\n", "    --no-penalty",
  //         "Turn off score penalty for repeat regions (def=penalty on).");
  fprintf(out_, "  %-30s %s\n", "    --blast-path <path>",
          "Path to directory with blastpgp executable (or set BLAST_PATH)");
}

int CSBlastApp::Run() {
  int status = 0;
  Init();

  for (SeqVec::iterator it = queries_.begin(); it != queries_.end(); ++it) {
    PrepareForRun(*it);
    CSBlastIteration itr(opts_.iterations);

    while (itr) {
      LOG(INFO) << strprintf("Starting iteration %i ...", itr.IterationNumber());

      SavePssm();

      // Set number of output descriptions and alignments
      if (itr.IterationNumber() == opts_.iterations) {
        opts_.csblast['v'] = strprintf("%i", opts_.ndescr);
        opts_.csblast['b'] = strprintf("%i", opts_.nalis);
        csblast_->set_options(opts_.csblast);
      }

      // Run one iteration of CS-BLAST
      FILE* fout = opts_.outfile.empty() ? out_ :
        fopen(opts_.outfile.c_str(), it != queries_.begin() ? "a" : "w");
      if (!fout) throw Exception("Unable to write to '%s'!", opts_.outfile.c_str());
      BlastHits hits;
      status = csblast_->Run(fout, &hits);
      if (!opts_.outfile.empty()) fclose(fout);

      // Don't bother parsing the results if this is the last iteration
      if (status != 0 || opts_.iterations == 1 || opts_.emulate || hits.empty()) break;

      hits.Filter(opts_.inclusion);
      if (!hits.empty() && !hits[0].hsps.empty())
        ali_->Merge(Alignment<AA>(hits, opts_.best));
      LOG(INFO) << strprintf("Found %zu seqs in iteration %i (E-value < %5.0E)",
                             hits.size(), itr.IterationNumber(), opts_.inclusion);
      itr.Advance(hits);

      if (itr) {
        CountProfile<AA> ali_profile(*ali_, !opts_.global_weights);
        if (opts_.pc_neff == 0.0) {
          CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
          pssm_.reset(new Pssm(*it, pc_->AddTo(ali_profile, admix)));
        } else {
          pssm_.reset(new Pssm(*it, pc_->AddTo(ali_profile, opts_.pc_neff)));
        }
        csblast_->set_pssm(pssm_.get());
      }
    }

    // Save alignment of hits if results were in -m0 format
    if (opts_.csblast.find('R') == opts_.csblast.end() ||
        opts_.csblast['m'] == "0")
      SaveAlignment();
  }

  return status;
}

void CSBlastApp::Init() {
  // Read query sequences
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
  ReadAll(fin, queries_);
  fclose(fin);
  if (queries_.empty())
    throw Exception("No sequences found in '%s'!", opts_.infile.c_str());

  // Remove all but first sequence if starting from checkpoint or alignment
  if (opts_.csblast.find('R') != opts_.csblast.end() ||
      !opts_.ali_infile.empty() ||
      opts_.iterations > 1) {
    queries_.erase(queries_.begin() + 1, queries_.end());
  }

  // Setup pseudocount engine
  if (opts_.pc_engine == "lib") {
    fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    lib_.reset(new ContextLibrary<AA>(fin));
    fclose(fin);

    TransformToLog(*lib_);
    pc_.reset(new LibraryPseudocounts<AA>(*lib_,
                                          opts_.weight_center,
                                          opts_.weight_decay));

  } else if (opts_.pc_engine == "crf") {
    fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    crf_.reset(new Crf<AA>(fin));
    fclose(fin);

    pc_.reset(new CrfPseudocounts<AA>(*crf_));
  } else if (opts_.pc_engine == "lr") {
    fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    lrparams_.reset(new LrParams<AA>(fin));
    fclose(fin);

    pc_.reset(new LrPseudocounts<AA>(*lrparams_));
  } else {
    throw Exception("Unknown pseudocount engine '%s'!", opts_.pc_engine.c_str());
  }

  // if (!opts_.no_penalty) {
  //   // Setup repeat penalizer
  //   penalizer_.reset(new RepeatPenalizer<AA>(opts_.penalty_alpha,
  //                                                   opts_.penalty_beta,
  //                                                   opts_.penalty_score_min));
  //   pc_->set_repeat_penalizer(penalizer_.get());
  // }
}

void CSBlastApp::PrepareForRun(const Sequence<AA>& query) {
  // Setup PSSM of query profile with context-specific pseudocounts if no
  // restart file is provided
  if (opts_.csblast.find('R') == opts_.csblast.end()) {
    if (opts_.ali_infile.empty()) {
      if (opts_.pc_neff == 0.0) {
        ConstantAdmix admix(opts_.pc_admix);
        pssm_.reset(new Pssm(query, pc_->AddTo(query, admix)));
      } else {
        pssm_.reset(new Pssm(query, pc_->AddTo(query, opts_.pc_neff)));
      }

    } else {
      FILE* fp = fopen(opts_.ali_infile.c_str(), "r");
      Alignment<AA> query_ali(fp, PSI_ALIGNMENT);
      query_ali.AssignMatchColumnsBySequence(0);
      fclose(fp);
      CountProfile<AA> ali_profile(query_ali, !opts_.global_weights);
      if (opts_.pc_neff == 0.0) {
        CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
        pssm_.reset(new Pssm(query, pc_->AddTo(ali_profile, admix)));
      } else {
        pssm_.reset(new Pssm(query, pc_->AddTo(ali_profile, opts_.pc_neff)));
      }
    }
  }

  // Reset number of output hits and alignments
  opts_.csblast['v']  = strprintf("%i", kNumOutputAlis);
  opts_.csblast['b']  = strprintf("%i", kNumOutputAlis);

  // Setup CS-BLAST engine
  if (opts_.csblast.find('R') == opts_.csblast.end())
    csblast_.reset(new CSBlast(&query, pssm_.get(), opts_.csblast));
  else
    csblast_.reset(new CSBlast(&query, opts_.csblast));

  // Set path to PSI-BLAST executable
  if (!opts_.blast_path.empty())
    csblast_->set_exec_path(opts_.blast_path);

  // Set BLAST call emulation
  csblast_->set_emulate(opts_.emulate);

  // Setup alignment of included sequences
  ali_.reset(new Alignment<AA>(query));
}

void CSBlastApp::SavePssm() const {
  if (!opts_.checkpointfile.empty() && pssm_) {
    FILE* fchk = fopen(opts_.checkpointfile.c_str(), "wb");
    if (!fchk)
      throw Exception("Can't write checkpoint '%s'!", opts_.checkpointfile.c_str());
    pssm_->Write(fchk);
    fclose(fchk);
  }
}

void CSBlastApp::SaveAlignment() const {
  if (!opts_.ali_outfile.empty()) {
    FILE* fali = fopen(opts_.ali_outfile.c_str(), "w");
    if (!fali) throw Exception("Can't write to '%s'!", opts_.ali_outfile.c_str());
    ali_->Write(fali, PSI_ALIGNMENT);
    fclose(fali);
  }
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSBlastApp().main(argc, argv, stdout, "csblast");
}
