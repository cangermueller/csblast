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

struct CSTrainsetNeffAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSTrainsetNeffAppOptions() { Init(); }
  virtual ~CSTrainsetNeffAppOptions() {}

  // Set csbuild default parameters
  void Init() {
    n = -1;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input count profile.
  string infile;
  // The output file with Neff of the training profiles.
  string outfile_x;
  // The output file with Neff of the training profiles' central column.
  string outfile_y;
  // Maximum number of training samples
  int n;

};  // CSTrainsetNeffAppOptions


template<class Abc>
class CSTrainsetNeffApp : public Application {
 private:
  // Runs the count_profile_neff application.
  virtual int Run();
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& ops);
  // Prints options summary to stream.
  virtual void PrintOptions() const;
  // Prints short application description.
  virtual void PrintBanner() const;
  // Prints usage banner to stream.
  virtual void PrintUsage() const;
  // Writes Neff to output file.
  inline void Write(FILE* fout, double neff) const;

  // Parameter wrapper
  CSTrainsetNeffAppOptions opts_;
};  // class CSTrainsetNeffApp



template<class Abc>
void CSTrainsetNeffApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('x', "outfile-x", opts_.outfile_x, opts_.outfile_x);
  ops >> Option('y', "outfile-y", opts_.outfile_y, opts_.outfile_y);
  ops >> Option('N', "size", opts_.n, opts_.n);
  opts_.Validate();
}

template<class Abc>
void CSTrainsetNeffApp<Abc>::PrintBanner() const {
  fputs("Calculates Neff in training profiles of a training set.\n", out_);
}

template<class Abc>
void CSTrainsetNeffApp<Abc>::PrintUsage() const {
  fputs("Usage: cstrainset_neff -i <infile> [options]\n", out_);
}

template<class Abc>
void CSTrainsetNeffApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input training set");
  fprintf(out_, "  %-30s %s\n", "-x, --outfile-x <file>",
          "Output file with Neff in the training profiles");
  fprintf(out_, "  %-30s %s\n", "-y, --outfile-y <file>",
          "Output file with Neff in the column to be predicted");
  fprintf(out_, "  %-30s %s (def=%i)\n", "-N, --size",
          "Maximum number of training samples", opts_.n);
}

template<class Abc>
int CSTrainsetNeffApp<Abc>::Run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Can't read input file '%s'!", opts_.infile.c_str());

  FILE* fout_x = NULL;
  if (opts_.outfile_x.length()) {
    fout_x = fopen(opts_.outfile_x.c_str(), "w");
    if (!fout_x)
      throw Exception("Can't open output file '%s'!", opts_.outfile_x.c_str());
  }

  FILE* fout_y = NULL;
  if (opts_.outfile_y.length()) {
    fout_y = fopen(opts_.outfile_y.c_str(), "w");
    if (!fout_y)
      throw Exception("Can't open output file '%s'!", opts_.outfile_y.c_str());
  }

  double neff_x = 0.0;
  double neff_y = 0.0;
  fprintf(out_, "Computing Neff...\n");
  if (GetFileExt(opts_.infile) == "tsq") {
    vector<TrainingSequence<Abc> > trainset;
    ReadAll(fin, trainset, opts_.n);
    for (size_t i = 0; i < trainset.size(); ++i) {
      double nx = 1.0;
      double ny = 0.0;
      for (size_t a = 0; a < Abc::kSize; ++a)
        ny += trainset[i].y[a];
      neff_y += ny;
      Write(fout_x, nx);
      Write(fout_y, ny);
    }
    neff_x = 1.0;
    neff_y /= trainset.size();

  } else {
    vector<TrainingProfile<Abc> > trainset;
    ReadAll(fin, trainset, opts_.n);
    for (size_t i = 0; i < trainset.size(); ++i) {
      double nx = Neff(trainset[i].x);
      double ny = 0;
      for (size_t a = 0; a < Abc::kSize; ++a)
        ny += trainset[i].y[a];
      neff_x += nx;
      neff_y += ny;
      Write(fout_x, nx);
      Write(fout_y, ny);
    }
    neff_x /= trainset.size();
    neff_y /= trainset.size();
  }

  fclose(fin);
  if (fout_x) fclose(fout_x);
  if (fout_y) fclose(fout_y);
  fprintf(out_, "Neff_x = %1.5f\n", neff_x);
  fprintf(out_, "Neff_y = %1.5f\n", neff_y);
  return 0;
}

template<class Abc>
inline void CSTrainsetNeffApp<Abc>::Write(FILE* fout, double neff) const {
  if (fout) fprintf(fout, "%8.5f\n", neff);
}



}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSTrainsetNeffApp<cs::Dna>().main(argc, argv, stdout, "cstrainset_neff");
  else
    return cs::CSTrainsetNeffApp<cs::AA>().main(argc, argv, stdout, "cstrainset_neff");
}
