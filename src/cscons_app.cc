// Copyright 2012, Christof Angermueller

#include "cs.h"
#include "alignment-inl.h"
#include "application.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

template<class Abc>
class ConsScoresEngine {
  public:
    virtual ~ConsScoresEngine() {};
    virtual vector<double> operator() (const Profile<Abc>& profile) const = 0;
};

class BasicConsScoresEngine : public ConsScoresEngine<AA> {
  public:
    BasicConsScoresEngine() {
      ss_.reset(new BlosumMatrix(BLOSUM62));
    }

    vector<double> operator() (const Profile<AA>& profile) const {
      vector<double> scores(profile.length());
      for (size_t i = 0; i < profile.length(); ++i) {
        double nom = 0.0;
        double denom = 0.0;
        for (size_t a = 0; a < AA::kSize; ++a) {
          nom += SQR(profile[i][a]) / ss_->p(a);
          denom += profile[i][a] / ss_->p(a);
        }
        scores[i] = log(nom) / log(denom);
        assert(scores[i] >= 0.0 && scores[i] <= 1.0);
      }
      return scores;
    }

  private:
    scoped_ptr<SubstitutionMatrix<AA> > ss_;
};

struct CSConsAppOptions {
  CSConsAppOptions() { Init(); }

  void Init() {
  }

  // Validates parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }
  
  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
};  // struct CSConsAppOptions


template <class Abc>
class CSConsApp : public Application {
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

  // Parameter wrapper
  CSConsAppOptions opts_;
};  // class CSConsApp



template <class Abc>
void CSConsApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);

  opts_.Validate();
}

template <class Abc>
void CSConsApp<Abc>::PrintBanner() const {
  fputs("Computes conservation scores for alignment columns.\n", out_);
}

template <class Abc>
void CSConsApp<Abc>::PrintUsage() const {
  fputs("Usage: cscons -i <infile> [options]\n", out_);
}

template <class Abc>
void CSConsApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with alignment");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with column specific conservation scores (def=stdout)");
}

template <class Abc>
int CSConsApp<Abc>::Run() {
  // Read input profile
  fprintf(out_, "Reading input profile from '%s'...\n", opts_.infile.c_str());
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read from '%s'!", opts_.infile.c_str());
  CountProfile<Abc> count_profile(fin);
  Profile<Abc> profile = count_profile.counts;
  Normalize(profile, 1.0);

  // Setup conservation scores engine
  scoped_ptr<ConsScoresEngine<Abc> >  cons_scores_engine;
  cons_scores_engine.reset(new BasicConsScoresEngine());

  // Computes conservation scores
  fprintf(out_, "Computing conservation scores for profile with %zu columns ...\n", profile.length());
  vector<double> cons_scores = (*cons_scores_engine)(profile);
  double mean = 0.0;
  for (size_t i = 0; i < cons_scores.size(); ++i) {
    mean += cons_scores[i];
  }
  mean /= cons_scores.size();
  double sd = 0.0;
  for (size_t i = 0; i < cons_scores.size(); ++i) {
    sd += SQR(cons_scores[i] - mean);
  }
  sd = sqrt(sd / (cons_scores.size() > 1 ? cons_scores.size() - 1 : cons_scores.size()));

  // Output conservation scores
  FILE* fout;
  if (opts_.outfile.empty()) {
    fout = out_;
  } else {
    fout = fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Unable to write to '%s'!", opts_.outfile.c_str());
    fprintf(out_, "Writing conservation scores to '%s' ...\n", opts_.outfile.c_str());
  }
  fprintf(fout, "# %6s  %8s  %8s\n", "NR", "SCORE", "Z-SCORE");
  for (size_t i = 0; i < cons_scores.size(); ++i) {
    fprintf(fout, "%8zu  %8.4f  %8.4f\n", i + 1, cons_scores[i], (cons_scores[i] - mean) / sd);
  }
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSConsApp<cs::AA>().main(argc, argv, stdout, "cscons");
}
