// Copyright 2011, Christof Angermueller

#include "cs.h"
#include "application.h"
#include "count_profile-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CountProfileNeffAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CountProfileNeffAppOptions() { Init(); }
  virtual ~CountProfileNeffAppOptions() {}

  // Set csbuild default parameters
  void Init() { }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input count profile.
  string infile;
};  // CountProfileNeffAppOptions


template<class Abc>
class CountProfileNeffApp : public Application {
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
  // Return the Neff in the count profile.
  inline double GetNeff(const CountProfile<Abc>& cp);

  // Parameter wrapper
  CountProfileNeffAppOptions opts_;
};  // class CountProfileNeffApp



template<class Abc>
void CountProfileNeffApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  opts_.Validate();
}

template<class Abc>
void CountProfileNeffApp<Abc>::PrintBanner() const {
  fputs("Calculate Neff of a count profile.\n", out_);
}

template<class Abc>
void CountProfileNeffApp<Abc>::PrintUsage() const {
  fputs("Usage: count_profile_neff -i <infile> [options]\n", out_);
}

template<class Abc>
void CountProfileNeffApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with count profile");
}

// Returns the Neff in given count profile.
template<class Abc>
inline double CountProfileNeffApp<Abc>::GetNeff(const CountProfile<Abc>& cp) {
  double neff_sum = 0.0;
  for (size_t i = 0; i < cp.length(); ++i) {
    for (size_t a = 0; a < Abc::kSize; ++a) {
      neff_sum += cp.counts[i][a] * log(cp.counts[i][a]);
    }
  }
  return exp(-neff_sum / cp.length());
}

template<class Abc>
int CountProfileNeffApp<Abc>::Run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Can't read input file '%s'!", opts_.infile.c_str());
  CountProfile<Abc> cp(fin);
  fclose(fin);
  double neff = GetNeff(cp);
  fprintf(out_, "Neff = %1.5f\n", neff);
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CountProfileNeffApp<cs::Dna>().main(argc, argv, stdout, "count_profile_neff");
  else
    return cs::CountProfileNeffApp<cs::AA>().main(argc, argv, stdout, "count_profile_neff");
}
