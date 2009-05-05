// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <string>
#include <vector>

#include "globals.h"
#include "alignment-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "getopt_pp.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct CSSampleAppOptions {
  CSSampleAppOptions()
      : sample_size(kMaxInt),
        window_length(0),
        sample_rate(0.2f) {}

  virtual ~CSSampleAppOptions() {}

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (outfile.empty()) throw Exception("No output file provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // The number of profiles to be sampled out.
  int sample_size;
  // The number of columns in the context window.
  int window_length;
  // Fraction of profile windows sampled from each full length alignment
  // or sequence.
  float sample_rate;
};  // CSSampleAppOptions


template<class Alphabet>
class CSSampleApp : public Application {
 private:
  typedef vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;
  typedef typename vector<int>::const_iterator index_iterator;

  // Runs the csbuild application.
  virtual int Run();
  // Parses command line options.
  virtual void parse_options(GetOpt_pp* options);
  // Prints options summary to stream.
  virtual void print_options() const;
  // Prints short application description.
  virtual void print_description() const;
  // Prints usage banner to stream.
  virtual void print_banner() const;
  // Samples profiles from database of profiles.
  void Sample();

  // Parameter wrapper
  CSSampleAppOptions opts_;
  // Database of profiles to sample from.
  profile_vector database_;
  // Sampled profiles.
  profile_vector samples_;
};  // CSSampleApp



template<class Alphabet>
void CSSampleApp<Alphabet>::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('N', "sample-size", opts_.sample_size,
                     opts_.sample_size);
  *options >> Option('W', "window-length", opts_.window_length,
                     opts_.window_length);
  *options >> Option('s', "sample-rate", opts_.sample_rate,
                     opts_.sample_rate);

  opts_.Validate();

  if (opts_.outfile.empty())
    opts_.outfile = get_file_basename(opts_.infile, false) + "prf";
}

template<class Alphabet>
void CSSampleApp<Alphabet>::print_description() const {
  fputs("Sample (context) profiles from a large profile database.\n",
        stream());
}

template<class Alphabet>
void CSSampleApp<Alphabet>::print_banner() const {
  fputs("Usage: cssample -i <infile> -o <outfile> [options]\n", stream());
}

template<class Alphabet>
void CSSampleApp<Alphabet>::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <file>",
          "Path to input file with profile database");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <file>",
          "Path for output file with sampled profiles");
  fprintf(stream(), "  %-30s %s\n", "-N, --num-profiles [0,inf[",
          "Maximal number of profiles to sample (def=inf)");
  fprintf(stream(), "  %-30s %s\n", "-W, --window-length [0,inf[",
          "Sample context profiles of length W instead of full-length profiles");
  fprintf(stream(), "  %-30s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]",
          "Fraction of context profiles sampled per full-length profile",
          opts_.sample_rate);
}

template<class Alphabet>
void CSSampleApp<Alphabet>::Sample() {
  fprintf(stream(), "Sampling %i profiles from pool of %i profiles ...\n",
          opts_.sample_size, static_cast<int>(database_.size()));
  fflush(stream());

  // Iterate over input data and build counts profiles either by full-length
  // conversion or by sampling of context windows
  for (profile_iterator it = database_.begin(); it != database_.end()
         && static_cast<int>(samples_.size()) < opts_.sample_size; ++it) {
    if (opts_.window_length == 0) {  // add full length profile to samples
      samples_.push_back(*it);

    } else if ((*it)->num_cols() >= opts_.window_length) {
      // Sample context windows if profile has sufficient length
      vector<int> idx;  // sample of indices
      for (int j = 0; j <= (*it)->num_cols() - opts_.window_length; ++j)
        idx.push_back(j);
      random_shuffle(idx.begin(), idx.end());
      const int sample_size = iround(opts_.sample_rate * idx.size());
      // sample only a fraction of the profile indices.
      idx.erase(idx.begin() + sample_size, idx.end());

      // Add sub-profiles at sampled indices to HMM
      for (index_iterator ii = idx.begin(); ii != idx.end()
             && static_cast<int>(samples_.size()) < opts_.sample_size; ++ii) {
        shared_ptr< CountProfile<Alphabet> > p(
            new CountProfile<Alphabet>(**it, *ii, opts_.window_length));
        samples_.push_back(p);
      }
    }
  }
}

template<class Alphabet>
int CSSampleApp<Alphabet>::Run() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Unable to read from input file '%s'!",
                    opts_.infile.c_str());
  fprintf(stream(), "Reading profiles from %s ...",
          get_file_basename(opts_.infile).c_str());
  fflush(stream());

  CountProfile<Alphabet>::readall(fin, &database_);

  fprintf(stream(), " %i profiles read\n", static_cast<int>(database_.size()));
  fclose(fin);

  // Sample profiles
  random_shuffle(database_.begin(), database_.end());
  Sample();

  // Write sampled profiles to outfile
  FILE* fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Unable to write profiles to output file '%s'!",
                    opts_.outfile.c_str());
  int num_cols = 0;
  for (profile_iterator it = samples_.begin(); it != samples_.end(); ++it) {
    (*it)->Write(fout);
    num_cols += (*it)->num_cols();
  }
  fprintf(stream(), "Wrote %i profiles with a total number of %i columns to %s\n",
          static_cast<int>(samples_.size()), num_cols, opts_.outfile.c_str());

  fclose(fout);

  return 0;
}

}  // namespace cs
