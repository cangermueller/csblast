// Copyright 2011, Christof Angermueller

#include "cs.h"
#include "application.h"
#include "count_profile-inl.h"
#include "training_sequence.h"
#include "training_profile.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct CSTrainProfilesAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSTrainProfilesAppOptions() { Init(); }
  virtual ~CSTrainProfilesAppOptions() {}

  // Set csbuild default parameters
  void Init() { }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (outfile.empty())
      outfile = GetBasename(infile, false) + ".tpr";
  }

  // The input file with training sequences.
  string infile;
  // The output file with training profiles.
  string outfile;

};  // CSTrainProfilesAppOptions


template<class Abc>
class CSTrainProfilesApp : public Application {
 private:
  // Runs the application.
  virtual int Run();
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& ops);
  // Prints options summary to stream.
  virtual void PrintOptions() const;
  // Prints short application description.
  virtual void PrintBanner() const;
  // Prints usage banner to stream.
  virtual void PrintUsage() const;

  // Parameter wrapper
  CSTrainProfilesAppOptions opts_;
};  // class CSTrainProfilesApp



template<class Abc>
void CSTrainProfilesApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  opts_.Validate();
}

template<class Abc>
void CSTrainProfilesApp<Abc>::PrintBanner() const {
  fputs("Converts a set of training sequences into a set of training profiles.\n", out_);
}

template<class Abc>
void CSTrainProfilesApp<Abc>::PrintUsage() const {
  fputs("Usage: cstrainprofiles -i <infile> [options]\n", out_);
}

template<class Abc>
void CSTrainProfilesApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with training sequences");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with training profiles (default: INFILE.tpr)");
}

template<class Abc>
int CSTrainProfilesApp<Abc>::Run() {
  size_t count = 0;
  FILE* fin;
  FILE* fout;
  fprintf(out_, "Converting training sequences to training profiles ...\n");
  fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
      throw Exception("Can't read training sequences from '%s'!", opts_.infile.c_str());
  fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
  while (!feof(fin)) {
    TrainingSequence<Abc> tseq(fin);
    TrainingProfile<Abc> tprof(tseq);
    tprof.Write(fout);
    count++;
    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
  fclose(fin);
  fclose(fout);
  fprintf(out_, "%zu training sequences converted to training profiles!\n", count);
  return 0;
}



}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSTrainProfilesApp<cs::Dna>().main(argc, argv, stdout, "cstrainprofiles");
  else
    return cs::CSTrainProfilesApp<cs::AA>().main(argc, argv, stdout, "cstrainprofiles");
}
