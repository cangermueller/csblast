/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cs.h"
#include "alignment-inl.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf_pseudocounts-inl.h"
#include "crf-inl.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "pssm.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSBuildAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSBuildAppOptions() { Init(); }
  virtual ~CSBuildAppOptions() {}

  // Set csbuild default parameters
  void Init() {
    informat         = "auto";
    outformat        = "prf";
    pc_admix         = 0.90;
    pc_neff          = 0.0;
    pc_ali           = 12.0;
    pc_engine        = "auto";
    match_assign     = kAssignMatchColsByQuery;
    global_weights   = false;
    weight_center    = 1.6;
    weight_decay     = 0.85;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (outformat == "chk" && match_assign != kAssignMatchColsByQuery)
      throw Exception("Can't write checkpoint if input alignment has no query!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Input file with context profile library or HMM
  string modelfile;
  // Input file format
  string informat;
  // Output file format
  string outformat;
  // Overall pseudocount admixture
  double pc_admix;
  // Target Neff for pseudocounts admixture
  double pc_neff;
  // Constant in pseudocount calculation for alignments
  double pc_ali;
  // Pseudocount engine
  string pc_engine;
  // Match column assignment for FASTA alignments
  int match_assign;
  // Use global instead of position specific weights for sequence weighting.
  bool global_weights;
  // Weight of central column in multinomial emission
  double weight_center;
  // Exponential decay of window weights
  double weight_decay;
};  // CSBuildAppOptions


template<class Abc>
class CSBuildApp : public Application {
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
  // Writes profile to outfile
  void WriteProfile(const CountProfile<Abc>& profile) const;
  // Writes PSI-BLAST checkpoint file
  void WriteCheckpoint(const Sequence<Abc>& query,
                       const CountProfile<Abc>& cp) const;

  // Parameter wrapper
  CSBuildAppOptions opts_;
  // Profile library for context pseudocounts
  scoped_ptr<ContextLibrary<Abc> > lib_;
  // CRF for CRF context pseudocounts
  scoped_ptr< Crf<Abc> > crf_;
  // Pseudocount engine
  scoped_ptr< Pseudocounts<Abc> > pc_;
};  // class CSBuildApp



template<class Abc>
void CSBuildApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('I', "informat", opts_.informat, opts_.informat);
  ops >> Option('O', "outformat", opts_.outformat, opts_.outformat);
  ops >> Option('M', "match-assign", opts_.match_assign, opts_.match_assign);
  ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  ops >> Option('z', "pc-neff", opts_.pc_neff, opts_.pc_neff);
  ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
  ops >> Option('p', "pc-engine", opts_.pc_engine, opts_.pc_engine);
  ops >> OptionPresent(' ', "global-weights", opts_.global_weights);
  ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
  ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);

  opts_.Validate();

  if (opts_.outfile.empty())
    opts_.outfile = GetBasename(opts_.infile, false) + ".prf";
  if (opts_.informat == "auto")
    opts_.informat = GetFileExt(opts_.infile);
  if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
    opts_.pc_engine = GetFileExt(opts_.modelfile);
}

template<class Abc>
void CSBuildApp<Abc>::PrintBanner() const {
  fputs("Build a profile or PSSM from an alignment or sequence.\n", out_);
}

template<class Abc>
void CSBuildApp<Abc>::PrintUsage() const {
  fputs("Usage: csbuild -i <infile> [options]\n", out_);
  fputs("       csbuild -i <infile> -D <contextdata> [options]\n", out_);
}

template<class Abc>
void CSBuildApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with alignment or sequence");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file for generated profile (def: <infile>.prf)");
  fprintf(out_, "  %-30s %s (def=%s)\n", "-I, --informat seq|fas|...",
          "Input format: seq, fas, a2m, or a3m", opts_.informat.c_str());
  fprintf(out_, "  %-30s %s (def=%s)\n", "-O, --outformat prf|chk",
          "Outformat: profile or PSI-BLAST checkpoint", opts_.outformat.c_str());
  fprintf(out_, "  %-30s %s\n", "-M, --match-assign [0:100]",
          "Make all FASTA columns with less than X% gaps match columns");
  fprintf(out_, "  %-30s %s\n", "", "(def: make columns with residue in "
          "first sequence match columns)");
  fprintf(out_, "  %-30s %s (def=off)\n", "-D, --context-data <file>",
          "Add context-specific pseudocounts with profile library");
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admixture for context-specific pseudocounts",
          opts_.pc_admix);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-z, --pc-neff [1,inf[",
          "Target Neff for pseudocounts admixture", opts_.pc_neff);
  fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant in pseudocount calculation for alignments", opts_.pc_ali);
  fprintf(out_, "  %-30s %s (def=%s)\n", "    --pc-engine crf|lib",
           "Specify engine for pseudocount generation", opts_.pc_engine.c_str());
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
         "Weight of central profile column for context-specific pseudocounts", opts_.weight_center);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
          "Parameter for exponential decay of window weights", opts_.weight_decay);
  fprintf(out_, "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific sequence weights (def=off)");
}

template<class Abc>
void CSBuildApp<Abc>::WriteProfile(const CountProfile<Abc>& profile) const {
  FILE* fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Can't write to output file '%s'!", opts_.outfile.c_str());
  profile.Write(fout);
  fprintf(out_, "Wrote profile with %zu columns to %s\n", profile.counts.length(),
          opts_.outfile.c_str());
  fclose(fout);
}

template<class Abc>
void CSBuildApp<Abc>::WriteCheckpoint(const Sequence<Abc>& /* query */,
                                      const CountProfile<Abc>& /* cp */) const {}

template<>
void CSBuildApp<AA>::WriteCheckpoint(const Sequence<AA>& query,
                                     const CountProfile<AA>& cp) const {
  Pssm pssm(query, cp.counts);
  Normalize(pssm.profile, 1.0);

  FILE* fout = fopen(opts_.outfile.c_str(), "wb");
  if (!fout) throw Exception("Unable to write '%s'!", opts_.outfile.c_str());
  pssm.Write(fout);
  fprintf(out_, "Wrote profile with %zu columns as checkpoint file to %s\n",
          cp.counts.length(), opts_.outfile.c_str());
  fclose(fout);
}

template<class Abc>
int CSBuildApp<Abc>::Run() {
  // Setup pseudocount engine
  if (!opts_.modelfile.empty() && opts_.pc_engine == "lib") {
    fprintf(out_, "Reading context library from %s ...\n",
            GetBasename(opts_.modelfile).c_str());
    FILE* fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    lib_.reset(new ContextLibrary<Abc>(fin));
    TransformToLog(*lib_);
    fclose(fin);
    pc_.reset(new LibraryPseudocounts<Abc>(*lib_, opts_.weight_center,
                                           opts_.weight_decay));
    pc_->SetTargetNeff(opts_.pc_neff);

  } else if (!opts_.modelfile.empty() && opts_.pc_engine == "crf") {
    fprintf(out_, "Reading CRF from %s ...\n",
            GetBasename(opts_.modelfile).c_str());
    FILE* fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    crf_.reset(new Crf<Abc>(fin));
    fclose(fin);
    pc_.reset(new CrfPseudocounts<Abc>(*crf_));
    pc_->SetTargetNeff(opts_.pc_neff);
  }

  CountProfile<Abc> profile;  // output profile
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Unable to read input file '%s'!", opts_.infile.c_str());

  if (opts_.informat == "seq") {  // build profile from sequence
    Sequence<Abc> seq(fin);
    profile = CountProfile<Abc>(seq);
    profile.name = GetBasename(opts_.infile, false);
    profile.name = profile.name.substr(0, profile.name.length() - 1);

    if (pc_) {
      fputs("Adding cs-pseudocounts ...\n", out_);
      ConstantAdmix admix(opts_.pc_admix);
      profile.counts = pc_->AddTo(seq, admix);
    }
    fprintf(out_, "Effective number of sequences exp(entropy) = %.2f\n", 
        Neff(profile.counts));

    if (opts_.outformat == "chk")
      WriteCheckpoint(seq, profile);
    else
      WriteProfile(profile);

  } else {  // build profile from alignment
    AlignmentFormat f = AlignmentFormatFromString(opts_.informat);
    Alignment<Abc> ali(fin, f);

    if (f == FASTA_ALIGNMENT) {
      if (opts_.match_assign == CSBuildAppOptions::kAssignMatchColsByQuery)
        ali.AssignMatchColumnsBySequence(0);
      else
        ali.AssignMatchColumnsByGapRule(opts_.match_assign);
    }
    profile = CountProfile<Abc>(ali, !opts_.global_weights);
    profile.name = GetBasename(opts_.infile, false);
    profile.name = profile.name.substr(0, profile.name.length() - 1);

    if (pc_) {
      fputs("Adding cs-pseudocounts ...\n", out_);
      CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
      profile.counts = pc_->AddTo(profile, admix);
      Normalize(profile.counts, profile.neff);
    }
    Profile<Abc> prof = profile.counts;
    Normalize(prof, 1.0);
    fprintf(out_, "Effective number of sequences exp(entropy) = %.2f\n", 
        Neff(prof));

    if (opts_.outformat == "chk")
      WriteCheckpoint(ali.GetSequence(0), profile);
    else
      WriteProfile(profile);
  }

  fclose(fin);
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSBuildApp<cs::Dna>().main(argc, argv, stdout, "csbuild");
  else
    return cs::CSBuildApp<cs::AA>().main(argc, argv, stdout, "csbuild");
}
