// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <memory>
#include <string>
#include <vector>

#include "globals.h"
#include "alignment-inl.h"
#include "amino_acid.h"
#include "blosum_matrix.h"
#include "clustering-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "getopt_pp.h"
#include "profile_library-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "nucleotide_matrix.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "substitution_matrix.h"
#include "utils-inl.h"

using namespace GetOpt;
using std::auto_ptr;
using std::string;
using std::vector;

namespace cs {

struct CSClustParams : public ClusteringParams {
  static const int kMatchColAssignByQuery = -1;

  CSClustParams()
      : num_profiles(0),
        lib_pseudocounts(1.0f),
        data_pseudocounts(0.01f),
        blosum_type("BLOSUM62"),
        nucleotide_match(1.0f),
        nucleotide_mismatch(-3.0f) { }

  virtual ~CSClustParams() { }

  // Validates the parameter settings and throws exception if needed.
  void validate() {
    if (num_profiles == 0 && libfile.empty())
      throw Exception("No value for number of profiles provided!");
    if (infile.empty())
      throw Exception("No input file with training data provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Directory for output and temporary files
  string directory;
  // Library input file for restarting
  string libfile;
  // The number of profiles in the profile library.
  int num_profiles;
  // Pseudocounts to be added to each library profile.
  float lib_pseudocounts;
  // Pseudocounts to be added to observed data counts.
  float data_pseudocounts;
  // BLOSUM matrix for pseudocount generation.
  string blosum_type;
  // Reward for a nucleotide match.
  float nucleotide_match;
  // Penalty for a nucleotide mismatch
  float nucleotide_mismatch;
};  // CSClustParams


template<class Alphabet>
class CSClustApp : public Application {
 private:
  typedef vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::iterator profile_iterator;

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
  // Prints substitution matrix options.
  void print_substitution_matrix_options() const;
  // Reads training data from infile.
  void read_training_data();
  // Initializes HMM from jumpstart file or by seeding from training data.
  void init_profile_library();
  // Initializes substitution matrix (specialized by alphabet type).
  void init_substitution_matrix();

  // Parameter wrapper
  CSClustParams params_;
  // Count profiles for training
  profile_vector data_;
  // Count profiles for training
  auto_ptr< ProfileLibrary<Alphabet> > lib_;
  // Substitution matrix for pseudocount generation
  auto_ptr< SubstitutionMatrix<Alphabet> > subst_matrix_;
};  // CSClustApp



template<class Alphabet>
void CSClustApp<Alphabet>::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", params_.infile, params_.infile);
  *options >> Option('o', "outfile", params_.outfile, params_.outfile);
  *options >> Option('d', "directory", params_.directory, params_.directory);
  *options >> Option('K', "num-profile", params_.num_profiles,
                     params_.num_profiles);
  *options >> Option('l', "likelihod-change", params_.log_likelihood_change,
                    params_.log_likelihood_change);
  *options >> Option('j', "jumpstart", params_.libfile, params_.libfile);
  *options >> Option('B', "blocks", params_.num_blocks, params_.num_blocks);
  *options >> Option('m', "matrix", params_.blosum_type, params_.blosum_type);
  *options >> Option('q', "mismatch-score", params_.nucleotide_mismatch,
                    params_.nucleotide_mismatch);
  *options >> Option('r', "match-score", params_.nucleotide_match,
                     params_.nucleotide_match);
  *options >> Option(' ', "data-pc", params_.data_pseudocounts,
                     params_.data_pseudocounts);
  *options >> Option(' ', "lib-pc", params_.lib_pseudocounts,
                     params_.lib_pseudocounts);
  *options >> Option(' ', "min-scans", params_.min_scans, params_.min_scans);
  *options >> Option(' ', "max-scans", params_.max_scans, params_.max_scans);
  *options >> Option(' ', "weight-center", params_.weight_center,
                     params_.weight_center);
  *options >> Option(' ', "weight-decay", params_.weight_decay,
                     params_.weight_decay);
  *options >> Option(' ', "epsilon", params_.epsilon_null,
                     params_.epsilon_null);
  *options >> Option(' ', "beta", params_.beta, params_.beta);

  params_.validate();

  if (!params_.directory.empty() && *params_.directory.rbegin() != kDirSep)
    params_.directory.append(1, kDirSep);
  if (params_.outfile.empty())
    params_.outfile = params_.directory +
      get_file_basename(params_.infile, false) + "lib";
}

template<class Alphabet>
void CSClustApp<Alphabet>::print_description() const {
  fputs("Cluster a training set of profile windows into a profile library.\n",
        stream());
}

template<class Alphabet>
void CSClustApp<Alphabet>::print_banner() const {
  fputs("Usage: csclust -i <infile> -K <num_profiles> [options]\n", stream());
  fputs("       csclust -i <infile> -j <libfile> [options]\n", stream());
}

template<class Alphabet>
void CSClustApp<Alphabet>::print_substitution_matrix_options() const {
  fprintf(stream(), "  %-30s %s (def=%.0f)\n", "-q, --mismatch-score <int>",
          "Penalty for a nucleotide mismatch", params_.nucleotide_mismatch);
  fprintf(stream(), "  %-30s %s (def=%.0f)\n",  "-r, --match-score <int>",
          "Reward for a nucleotide match", params_.nucleotide_match);
}

template<>
void CSClustApp<AminoAcid>::print_substitution_matrix_options() const {
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-m, --matrix <string>",
          "Substitution matrix: BLOSUM45, BLOSUM62, or BLOSUM80",
          params_.blosum_type.c_str());
}

template<class Alphabet>
void CSClustApp<Alphabet>::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Path to input file with profile windows");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Path for output file with profile library");
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-d, --directory <directory>",
          "Directory for temporary and output files",
          params_.directory.empty() ? "." : params_.directory.c_str());
  fprintf(stream(), "  %-30s %s\n", "-K, --num-profiles [0,inf[",
          "Number of profiles in the library");
  fprintf(stream(), "  %-30s %s (def=%3.1g)\n", "-l, --likelihood [0,inf[",
          "Maximal likelihood change per column for convergence",
          params_.log_likelihood_change);
  fprintf(stream(), "  %-30s %s\n", "-j, --jumpstart <filename>",
          "Jumpstart the clustering with a profile library");
  fprintf(stream(), "  %-30s %s\n", "-B, --blocks [0,N]",
          "Number of blocks for online training (def: B=N^3/8)");

  print_substitution_matrix_options();

  fprintf(stream(), "  %-30s %s (def=%i)\n", "    --min-scans [0,inf[",
          "Minimal number of training data scans", params_.min_scans);
  fprintf(stream(), "  %-30s %s (def=%i)\n", "    --max-scans [0,inf[",
          "Maximal number of training data scans", params_.max_scans);
  fprintf(stream(), "  %-30s %s (def=%3.1f)\n", "    --lib-pc [0,1]",
          "Pseudocounts for library profiles", params_.lib_pseudocounts);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --data-pc [0,1]",
          "Pseudocounts for training data", params_.data_pseudocounts);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --weight-center [0,1]",
          "Weight of central profile column in context window",
          params_.weight_center);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --weight-decay [0,1]",
          "Exponential decay of positional window weights",
          params_.weight_decay);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --epsilon [0,1]",
          "Start value for learning rate epsilon in online training",
          params_.epsilon_null);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --beta [0,1]",
          "Exponential decay of epsilon in online training", params_.beta);
}

template<class Alphabet>
void CSClustApp<Alphabet>::init_substitution_matrix() {
  subst_matrix_ =  auto_ptr< SubstitutionMatrix<Alphabet> >(
      new NucleotideMatrix(params_.nucleotide_match,
                           params_.nucleotide_mismatch));
}

template<>
void CSClustApp<AminoAcid>::init_substitution_matrix() {
  BlosumMatrix::Type type = blosum_matrix_type_from_string(params_.blosum_type);
  subst_matrix_ = auto_ptr< SubstitutionMatrix<AminoAcid> >(new BlosumMatrix(type));
}

template<class Alphabet>
void CSClustApp<Alphabet>::read_training_data() {
  FILE* fin = fopen(params_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Unable to read from input file '%s'!",
                    params_.infile.c_str());

  // Read data counts directly from serialized count profiles
  fprintf(stream(), "Reading training profiles from %s ...",
          get_file_basename(params_.infile).c_str());
  fflush(stream());
  CountProfile<Alphabet>::readall(fin, &data_);
  fprintf(stream(), " %i profiles read\n", data_.size());

  fclose(fin);
}

template<class Alphabet>
void CSClustApp<Alphabet>::init_profile_library() {
  if (params_.libfile.empty()) {  // seed the library from training data
    fprintf(stream(), "Initializing profile library by sampling %i context "
            "profiles from training profiles ...", params_.num_profiles);
    fflush(stream());
    MatrixPseudocounts<Alphabet> pc(subst_matrix_.get());
    SamplingProfileInitializer<Alphabet> profile_init(data_,
                                                      &pc,
                                                      params_.lib_pseudocounts);
    lib_ = auto_ptr< ProfileLibrary<Alphabet> >(
        new ProfileLibrary<Alphabet>(params_.num_profiles,
                                     data_[0]->num_cols(),
                                     profile_init));
    lib_->transform_to_logspace();
    fputc('\n', stream());

  } else {  // read profile library from jumpstart file
    FILE* fin = fopen(params_.libfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read from jumpstart file '%s'!",
                      params_.libfile.c_str());
    fprintf(stream(), "Reading profile library from %s ...",
            get_file_basename(params_.libfile).c_str());
    fflush(stream());
    lib_ = auto_ptr< ProfileLibrary<Alphabet> >(new ProfileLibrary<Alphabet>(fin));
    fputc('\n', stream());
    fclose(fin);
  }
}

template<class Alphabet>
int CSClustApp<Alphabet>::run() {
  init_substitution_matrix();
  read_training_data();
  init_profile_library();

  // Add pseudocounts to training data
  fprintf(stream(), "Adding pseudocounts to training profiles (admix=%.2f) ...",
          params_.data_pseudocounts);
  fflush(stream());
  MatrixPseudocounts<Alphabet> pc(subst_matrix_.get());
  for (profile_iterator ci = data_.begin(); ci != data_.end(); ++ci) {
    pc.add_to_profile(ConstantAdmixture(params_.data_pseudocounts), ci->get());
  }
  fputc('\n', stream());

  // Run EM clustering
  fprintf(stream(), "Clustering training data (K=%i, W=%i, N=%i) ...",
          lib_->num_profiles(), lib_->num_cols(), data_.size());
  fflush(stream());
  fputs("\n\n", stream());
  Clustering<Alphabet, CountProfile> clust(params_, data_, *lib_, stream());
  clust.run();

  // Write profile library to outfile
  FILE* fout = fopen(params_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Unable to write profiles to output file '%s'!",
                    params_.outfile.c_str());
  lib_->write(fout);
  fclose(fout);
  fprintf(stream(), "\nWrote profile library to %s\n", params_.outfile.c_str());

  return 0;
}

}  // namespace cs
