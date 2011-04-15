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
    informat = "auto";
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Input file with abstract state library
  string libfile;
  // Input file format
  string informat;

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

  opts_.Validate();

  if (opts_.informat == "auto")
    opts_.informat = GetFileExt(opts_.infile);
  if (opts_.outfile.empty())
    opts_.outfile = GetBasename(opts_.infile, false) + ".pdf";
  if (GetDirname(opts_.outfile).empty())
    opts_.outfile = "./" + opts_.outfile;
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
    profile_writer.WriteToFile(opts_.outfile);
    fprintf(out_, "Wrote output PDF to %s\n", opts_.outfile.c_str());

  } else if (opts_.informat == "ap62") {
    FILE* fin = fopen(opts_.infile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
    CountProfile<AS62> profile(fin);
    fclose(fin);

    fputs("Generating PDF with profile logos ...\n", out_);
    StateProfilePdfWriter<AS62, Abc> profile_writer(profile.counts, *lib_);
    profile_writer.WriteToFile(opts_.outfile);
    fprintf(out_, "Wrote output PDF to %s\n", opts_.outfile.c_str());
  }

  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSVizApp<cs::AA>().main(argc, argv, stdout, "csviz");
}
