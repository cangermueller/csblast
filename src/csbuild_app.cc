// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "alignment-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "getopt_pp.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSBuildParams {
  static const int kMatchColAssignByQuery = -1;

  CSBuildParams()
      : format("auto"),
        matchcol_assignment(kMatchColAssignByQuery),
        global_weights(false) { }

  virtual ~CSBuildParams() { }

  // Validates the parameter settings and throws exception if needed.
  void validate() {
    if (infile.empty()) throw Exception("No input file provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // File format of input alignment
  string format;
  // Match column assignment for FASTA alignments
  int matchcol_assignment;
  // Use global instead of position specific weights for profile construction.
  bool global_weights;
};  // CSBuildParams


template<class Alphabet>
class CSBuildApp : public Application {
 private:
  // Runs the csbuild application.
  virtual int run();
  // Parses command line options.
  virtual void parse_options(GetOpt_pp* options);
  // Prints options summary to stream.
  virtual void print_options() const;
  // Prints short application description.
  virtual void print_description() const;
  // Prints usage banner to stream.
  virtual void print_banner() const;

  CSBuildParams params_;  // parameter wrapper
};  // CSBuildApp



template<class Alphabet>
void CSBuildApp<Alphabet>::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", params_.infile, params_.infile);
  *options >> Option('o', "outfile", params_.outfile, params_.outfile);
  *options >> Option('f', "format", params_.format, params_.format);
  *options >> Option('M', "matchcol", params_.matchcol_assignment,
                     params_.matchcol_assignment);
  *options >> OptionPresent(' ', "global-weights", params_.global_weights);

  params_.validate();

  if (params_.outfile.empty())
    params_.outfile = get_file_basename(params_.infile, false) + "prf";
  if (params_.format == "auto")
    params_.format = get_file_ext(params_.infile);
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_description() const {
  fputs("Build a sequence profile from an alignment or sequence.\n",
        stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_banner() const {
  fputs("Usage: csbuild -i <infile> [options]\n", stream());
}

template<class Alphabet>
void CSBuildApp<Alphabet>::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Path to input file with alignment or sequence");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Path for output file with serialized profile");
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-f, --format <string>",
          "Input data format: seq, fas, a2m, or a3m", params_.format.c_str());
  fprintf(stream(), "  %-30s %s\n", "-M, --matchcol [0:100]",
          "Make all FASTA columns with less than X% gaps match columns");
  fprintf(stream(), "  %-30s %s\n", "", "(def: make columns with residue in "
          "first sequence match columns)");
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific weights for profiles");
}

template<class Alphabet>
int CSBuildApp<Alphabet>::run() {
  FILE* fin = fopen(params_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Unable to read from input file '%s'!",
                    params_.infile.c_str());
  FILE* fout = fopen(params_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Unable to write profiles to output file '%s'!",
                    params_.outfile.c_str());

  if (params_.format == "seq") {  // build profile from sequence
    Sequence<Alphabet> seq(fin);
    CountProfile<Alphabet> profile(seq);
    profile.write(fout);
    fprintf(stream(), "Wrote profile with %i columns to %s\n",
            profile.num_cols(), params_.outfile.c_str());

  } else {  // build profile from alignment
    typename Alignment<Alphabet>::Format f =
      alignment_format_from_string<Alphabet>(params_.format);
    Alignment<Alphabet> ali(fin, f);
    if (f == Alignment<Alphabet>::FASTA) {
      if (params_.matchcol_assignment == CSBuildParams::kMatchColAssignByQuery)
        ali.assign_match_columns_by_sequence(0);
      else
        ali.assign_match_columns_by_gap_rule(params_.matchcol_assignment);
    }
    CountProfile<Alphabet> profile(ali, !params_.global_weights);
    profile.write(fout);
    fprintf(stream(), "Wrote profile with %i columns to %s\n",
            profile.num_cols(), params_.outfile.c_str());
  }

  fclose(fout);
  fclose(fin);

  return 0;
}

}  // namespace cs
