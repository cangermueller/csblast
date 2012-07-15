// Copyright 2010, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "blosum_matrix.h"
#include "co_emission.h"
#include "context_library-inl.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "pdf_writer-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSVizAppOptions {
  static const int kAssignMatchColsByQuery = -1;

  CSVizAppOptions() { Init(); }
  virtual ~CSVizAppOptions() {}

  // Set csbuild default parameters
  void Init() {
    informat     = "auto";
    keep         = false;
    sort         = false;
    states       = "";
    max_states   = 0;
    crf_weights  = false;
    external_dir = "";
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  vector<size_t> States() const {
    vector<size_t> states_int;
    if (!states.empty()) {
      vector<std::string> s;
      Tokenize(states, ',', &s);
      if (s.empty()) {
        s.push_back(states);
      }
      for (size_t i = 0; i < s.size(); ++i) {
        vector<std::string> ss;
        Tokenize(s[i], "-", &ss);
        if (ss.size() == 2) {
          size_t first = atoi(ss[0].c_str());
          size_t last = atoi(ss[1].c_str());
          for (size_t i = first; i <= last; ++i) {
            states_int.push_back(i);
          }
        } else {
          states_int.push_back(atoi(s[i].c_str()));
        }
      }
    }
    return states_int;
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Input file with abstract state library
  string libfile;
  // Input file format
  string informat;
  // Keep temporary files
  bool keep;
  // Sort states
  bool sort;
  // Expression for selecting states
  string states;
  // Maximum number of states to be printed
  size_t max_states;
  // Visualize CRF weights instead of probabilities
  bool crf_weights;
  // Create figures in external directory
  string external_dir;

};  // CSVizAppOptions


template<class Abc>
class CSVizApp : public Application {
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
  // Prepares CRF states for visualization.
  std::vector<CrfState<Abc> > PrepareCrfStates(const Crf<Abc>& crf) const;

  // Parameter wrapper
  CSVizAppOptions opts_;
  // Profile library with abstract states
  scoped_ptr<ContextLibrary<Abc> > lib_;
};  // class CSVizApp



template<class Abc>
void CSVizApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('A', "alphabet", opts_.libfile, opts_.libfile);
  ops >> Option('I', "informat", opts_.informat, opts_.informat);
  ops >> OptionPresent('k', "keep", opts_.keep);
  ops >> OptionPresent(' ', "sort", opts_.sort);
  ops >> Option(' ', "states", opts_.states, opts_.states);
  ops >> Option(' ', "max-states", opts_.max_states, opts_.max_states);
  // ops >> OptionPresent(' ', "crf-weights", opts_.crf_weights);
  ops >> Option(' ', "external", opts_.external_dir, opts_.external_dir);

  opts_.Validate();

  if (opts_.informat == "auto")
    opts_.informat = GetFileExt(opts_.infile);
  if (opts_.outfile.empty())
    opts_.outfile = GetBasename(opts_.infile, false) + ".pdf";
  if (GetDirname(opts_.outfile).empty())
    opts_.outfile = "./" + opts_.outfile;
  if (*opts_.external_dir.rbegin() == kDirSep) 
    opts_.external_dir.substr(0, opts_.external_dir.size() - 1);
}

template<class Abc>
void CSVizApp<Abc>::PrintBanner() const {
  fputs("Draw profile logos for an amino acid or abstract state profile.\n", out_);
}

template<class Abc>
void CSVizApp<Abc>::PrintUsage() const {
  fputs("Usage: csviz -i <infile> [options]\n", out_);
}

template<class Abc>
void CSVizApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with amino acid or abstract state profile");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "PDF output file with profile logos (def: <infile>.pdf)");
  fprintf(out_, "  %-30s %s (def=%s)\n", "-I, --informat prf|ap62",
          "Input format: amino acid (prf) or abstract state profile (ap62)",
          opts_.informat.c_str());
  fprintf(out_, "  %-30s %s (def=off)\n", "-A, --alphabet <file>",
          "Abstract state alphabet consisting of exactly 62 states");
  fprintf(out_, "  %-30s %s (def=off)\n", "-k, --keep", "Keep temporary files");
  fprintf(out_, "  %-30s %s (def=off)\n", "    --sort", "Sort states by probability");
  fprintf(out_, "  %-30s %s (def=off)\n", "    --states <s[-s][,s[-s]]+>", 
          "Expression for selecting states");
  fprintf(out_, "  %-30s %s (def=all)\n", "    --max-states <int>", 
          "Maximum number of states to be printed");
  // fprintf(out_, "  %-30s %s (def=off)\n", "    --crf-weights", 
  //     "Visualize CRF weights instead of probabilities");
  fprintf(out_, "  %-30s %s\n", "    --external <directory>", 
      "Create figures in external directory");
}

template<class Abc>
int CSVizApp<Abc>::Run() {
  // Reading abstract state library
  if (!opts_.libfile.empty()) {
    fprintf(out_, "Reading abstract state alphabet from %s ...\n",
            GetBasename(opts_.libfile).c_str());
    FILE* fin = fopen(opts_.libfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.libfile.c_str());
    lib_.reset(new ContextLibrary<Abc>(fin));
    TransformToLin(*lib_);
    fclose(fin);
  }

  // Read input profile and generate output PDF
  if (opts_.informat == "prf") {
    FILE* fin = fopen(opts_.infile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
    CountProfile<Abc> profile(fin);
    fclose(fin);

    fputs("Generating PDF with profile logos ...\n", out_);
    ProfilePdfWriter<Abc> profile_writer(profile.counts);
    profile_writer.WriteToFile(opts_.outfile, opts_.keep);
    fprintf(out_, "Wrote output PDF to %s\n", opts_.outfile.c_str());

  } else if (opts_.informat == "ap62") {
    FILE* fin = fopen(opts_.infile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
    CountProfile<AS62> profile(fin);
    fclose(fin);

    fputs("Generating PDF with profile logos ...\n", out_);
    StateProfilePdfWriter<AS62, Abc> profile_writer(profile.counts, *lib_);
    profile_writer.WriteToFile(opts_.outfile, opts_.keep);
    fprintf(out_, "Wrote output PDF to %s\n", opts_.outfile.c_str());

  } else if (opts_.informat == "crf") {
    FILE* fin = fopen(opts_.infile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
    Crf<Abc> crf(fin);
    fclose(fin);

    fputs("Generating PDF from CRF ...\n", out_);
    vector<CrfState<Abc> > crf_states = PrepareCrfStates(crf);
    scoped_ptr<CrfStatePdfWriter<Abc> > crf_state_writer;
    if (opts_.crf_weights) {
      crf_state_writer.reset(new WeightCrfStatePdfWriter<Abc>());
    } else {
      crf_state_writer.reset(new ProbCrfStatePdfWriter<Abc>());
    }
    CrfPdfWriter<Abc> crf_writer(crf_states, *crf_state_writer);
    crf_writer.external_dir = opts_.external_dir;
    crf_writer.WriteToFile(opts_.outfile, opts_.keep);
    fprintf(out_, "Wrote output PDF to %s\n", opts_.outfile.c_str());
  }

  return 0;
}

template<class Abc>
std::vector<CrfState<Abc> > CSVizApp<Abc>::PrepareCrfStates(const Crf<Abc>& crf) const {
  std::vector<CrfState<Abc> > crf_states;

  std::vector<size_t> idx;
  if (opts_.states.empty()) {
    for (size_t i = 1; i <= crf.size(); ++i) {
      idx.push_back(i);
    }
  } else {
    idx = opts_.States();
  }
  for (size_t i = 0; i < idx.size(); ++i) {
    CrfState<Abc> crf_state = crf[idx[i] - 1];
    crf_state.name = strprintf("%zu: %+.1f", idx[i], crf_state.bias_weight);
    crf_states.push_back(crf_state);
  }

  if (opts_.sort) {
    std::sort(crf_states.begin(), crf_states.end());
    std::reverse(crf_states.begin(), crf_states.end());
  }

  if (opts_.max_states > 0 && crf_states.size() > opts_.max_states) {
    crf_states.erase(crf_states.begin() + opts_.max_states, crf_states.end());
  }

  return crf_states;
}


}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSVizApp<cs::AA>().main(argc, argv, stdout, "csviz");
}
