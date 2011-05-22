// Copyright 2011, Christof Angermueller

#include "cs.h"
#include "application.h"
#include "count_profile-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSCpNeffAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSCpNeffAppOptions() { Init(); }
  virtual ~CSCpNeffAppOptions() {}

  // Set csbuild default parameters
  void Init() { 
    neff_col = false;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input count profile.
  string infile;
  // Use the Neff per column for calculating the overall Neff
  bool neff_col;
};  // CSCpNeffAppOptions


template<class Abc>
class CSCpNeffApp : public Application {
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

  // Parameter wrapper
  CSCpNeffAppOptions opts_;
};  // class CSCpNeffApp



template<class Abc>
void CSCpNeffApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> OptionPresent('c', "neff-col", opts_.neff_col);
  opts_.Validate();
}

template<class Abc>
void CSCpNeffApp<Abc>::PrintBanner() const {
  fputs("Calculate Neff of a count profile.\n", out_);
}

template<class Abc>
void CSCpNeffApp<Abc>::PrintUsage() const {
  fputs("Usage: cscp_neff -i <infile> [options]\n", out_);
}

template<class Abc>
void CSCpNeffApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with count profile");
  fprintf(out_, "  %-30s %s(def=off)\n", "-c, --neff-col",
          "Use the Neff per column for calculating the overall Neff");
}

template<class Abc>
int CSCpNeffApp<Abc>::Run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Can't read input file '%s'!", opts_.infile.c_str());
  CountProfile<Abc> cp(fin);
  fclose(fin);
  double neff;
  if (opts_.neff_col) {
    neff = Neff(cp);
  } else {
    Normalize(cp.counts, 1.0);
    neff = Neff(cp.counts);
  }
  fprintf(out_, "Neff = %1.5f\n", neff);
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSCpNeffApp<cs::Dna>().main(argc, argv, stdout, "cscp_neff");
  else
    return cs::CSCpNeffApp<cs::AA>().main(argc, argv, stdout, "cscp_neff");
}
