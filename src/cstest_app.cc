// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf-inl.h"
#include "count_profile-inl.h"
#include "func.h"
#include "getopt_pp.h"
#include "matrix_pseudocounts-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"

using namespace GetOpt;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;

namespace cs {

struct CSTestAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSTestAppOptions() { Init(); }
  virtual ~CSTestAppOptions() {}

  // Set csbuild default parameters
  void Init() {
    pc_admix         = 0.9;
    weight_center    = 1.6;
    weight_decay     = 0.85;
    blosum_type      = "BLOSUM62";
    equalize         = false;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (testfile.empty()) throw Exception("No test set provided!");
    if (modelfiles.empty()) throw Exception("No model files provided!");
  }

  // Input file with training set.
  string testfile;
  // The input alignment file with training data.
  vector<string> modelfiles;
  // Output file for progress table.
  string outfile;
  // Overall pseudocount admixture
  double pc_admix;
  // Weight of central column in multinomial emission
  double weight_center;
  // Exponential decay of window weights
  double weight_decay;
  // BLOSUM matrix for pseudocount generation.
  string blosum_type;
  // Equalize all column counts to one?
  bool equalize;
};  // CSTestAppOptions


template<class Abc>
class CSTestApp : public Application {
 protected:
  typedef vector<TrainingSequence<Abc> > TrainingSet;
  typedef ContextLibFunc<Abc, TrainingSequence<Abc> > LibLikelihood;
  typedef CrfFunc<Abc, TrainingSequence<Abc> > CrfLikelihood;
  typedef pair<double, string> LoglikeNamePair;

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
  // Initializes substitution matrix.
  void InitSubstitutionMatrix();

  CSTestAppOptions opts_;
  TrainingSet testset_;
  scoped_ptr<SubstitutionMatrix<Abc> > sm_;
};  // class CSTestApp


template<class Abc>
void CSTestApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "testset", opts_.testfile, opts_.testfile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
  ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
  ops >> OptionPresent(' ', "equalize", opts_.equalize);
  ops >> Option(GetOpt_pp::EMPTY_OPTION, opts_.modelfiles);

  opts_.Validate();
}

template<class Abc>
void CSTestApp<Abc>::PrintBanner() const {
  fputs("Evaluate likelihood of context models on a test set.\n", out_);
}

template<class Abc>
void CSTestApp<Abc>::PrintUsage() const {
  fputs("Usage: cstest <modelfiles> -i <testset> [options]\n", out_);
}

template<class Abc>
void CSTestApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --testset <file>", "File with test set");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "Outfile with ranked models"
          " (def=stdout)");
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admixture for context library pseudocounts", opts_.pc_admix);
  fprintf(out_, "  %-30s %s (def=%4.2f)\n", "    --weight-center [0,inf[",
          "Weight of centralcolumn in context window", opts_.weight_center);
  fprintf(out_, "  %-30s %s (def=%4.2f)\n", "    --weight-decay [0,1]",
          "Exponential decay of positional window weights", opts_.weight_decay);
  fprintf(out_, "  %-30s %s\n", "    --equalize",
          "Equalize all column count sums to one (def=off)");
}

template<class Abc>
void CSTestApp<Abc>::InitSubstitutionMatrix() {
  BlosumType type = BlosumTypeFromString(opts_.blosum_type);
  sm_.reset(new BlosumMatrix(type));
}

template<class Abc>
int CSTestApp<Abc>::Run() {
  InitSubstitutionMatrix();

  // Read test set
  FILE* fin;
  fprintf(out_, "Reading test set from %s ...\n",
          GetBasename(opts_.testfile).c_str());
  fin = fopen(opts_.testfile.c_str(), "r");
  if (!fin) throw Exception("Can't read from '%s'!", opts_.testfile.c_str());
  ReadAll(fin, testset_);
  if (opts_.equalize) {
    for (size_t n = 0; n < testset_.size(); ++n)
      Normalize(&testset_[n].y[0], Abc::kSize);
  }
  fclose(fin);
  fprintf(out_, "%zu records read\n\n", testset_.size());

  // Test all model files and keep track of longest model name
  size_t max_len = 0;
  for (size_t m = 0; m < opts_.modelfiles.size(); ++m) {
    fin = fopen(opts_.modelfiles[m].c_str(), "r");
    if (!fin) throw Exception("Can't read from '%s'!", opts_.modelfiles[m].c_str());
    fclose(fin);
    max_len = MAX(GetBasename(opts_.modelfiles[m]).length(), max_len);
  }
  string model_str = string("Model") + string(max_len - 5, ' ');
  fprintf(out_, "%s  %-12s %8s\n", model_str.c_str(), "Progress", "Loglike");
  fprintf(out_, "%s\n", string(max_len + 23, '-').c_str());

  // Evaluate models on test set and keep track of likelihoods
  vector<LoglikeNamePair> ranking;
  for (size_t m = 0; m < opts_.modelfiles.size(); ++m) {
    string basename = GetBasename(opts_.modelfiles[m]);
    string tmp = basename + string(max_len - basename.length(), ' ');
    fprintf(out_, "%s  ", tmp.c_str());
    ProgressBar prog_bar(out_, 12);
    double loglike = 0.0;

    if (opts_.modelfiles[m].find(".lib") != std::string::npos) {
      fin = fopen(opts_.modelfiles[m].c_str(), "r");
      ContextLibrary<Abc> lib(fin);
      fclose(fin);
      TransformToLog(lib);

      LibLikelihood func(testset_, *sm_, opts_.weight_center,
                         opts_.weight_decay, opts_.pc_admix);
      prog_bar.Init(testset_.size() * lib.size());
      loglike = func(lib, &prog_bar);

    } else if (opts_.modelfiles[m].find(".crf") != std::string::npos) {
      fin = fopen(opts_.modelfiles[m].c_str(), "r");
      Crf<Abc> crf(fin);
      fclose(fin);

      CrfLikelihood func(testset_, *sm_);
      prog_bar.Init(testset_.size() * crf.size());
      loglike = func(crf, &prog_bar);
    }

    fprintf(out_, " %8.4f\n", loglike / testset_.size());
    ranking.push_back(make_pair(loglike / testset_.size(), basename));
  }
  fputs("\n", out_);

  // Rank models by likelihood
  sort(ranking.begin(), ranking.end());
  reverse(ranking.begin(), ranking.end());

  // Print models ranked by likelihood
  assert_eq(opts_.modelfiles.size(), ranking.size());
  FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");
  if (!fout) throw Exception("Can't write to file '%s'!", opts_.outfile.c_str());
  fprintf(fout, "%4s %s %8s\n", "Rank", model_str.c_str(), "Loglike");
  fprintf(fout, "%s\n", string(max_len + 14, '-').c_str());
  for (size_t m = 0; m < opts_.modelfiles.size(); ++m) {
    string basename = ranking[m].second;
    string tmp = basename + string(max_len - basename.length(), ' ');
    fprintf(fout, "%-4zu %s %8.4f\n", m + 1, tmp.c_str(), ranking[m].first);
  }
  if (!opts_.outfile.empty()) fclose(fout);

  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSTestApp<cs::AA>().main(argc, argv, stdout, "cstest");
}
