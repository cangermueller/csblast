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
#include "baum_welch_training-inl.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "getopt_pp.h"
#include "hmm-inl.h"
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

struct CSTrainAppOptions : public BaumWelchOptions {
  static const int kMatchColAssignByQuery = -1;

  CSTrainAppOptions()
      : format("auto"),
        matchcol_assignment(kMatchColAssignByQuery),
        num_states(0),
        window_length(1),
        sample_rate(0.2f),
        state_pseudocounts(1.0f),
        data_pseudocounts(0.01f),
        global_weights(false),
        blosum_type("BLOSUM62"),
        nucleotide_match(1.0f),
        nucleotide_mismatch(-3.0f) {}

  virtual ~CSTrainAppOptions() {}

  // Validates the parameter settings and throws exception if needed.
  void validate() {
    if (num_states == 0 && hmmfile.empty())
      throw Exception("No value for number of HMM states provided!");
    if (infile.empty())
      throw Exception("No input file with training data provided!");
  }

  // The input alignment file with training data.
  string infile;
  // The output file for the trained HMM.
  string outfile;
  // Directory for output and temporary files
  string directory;
  // File format of input alignment
  string format;
  // HMM input file for restarting
  string hmmfile;
  // Match column assignment for FASTA alignments
  int matchcol_assignment;
  // The number of states in the HMM to train.
  int num_states;
  // The number of columns in the context window.
  int window_length;
  // Fraction of profile windows sampled from each full length subject.
  float sample_rate;
  // Pseudocounts to be added to each state profile.
  float state_pseudocounts;
  // Pseudocounts to be added to observed data counts.
  float data_pseudocounts;
  // Use global instead of position specific weights for profile construction.
  bool global_weights;
  // BLOSUM matrix for pseudocount generation.
  string blosum_type;
  // Reward for a nucleotide match.
  float nucleotide_match;
  // Penalty for a nucleotide mismatch
  float nucleotide_mismatch;
};  // CSTrainAppOptions


template<class Alphabet>
class CSTrainApp : public Application {
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
  void init_hmm();
  // Initializes substitution matrix (specialized by alphabet type).
  void init_substitution_matrix();

  // Parameter wrapper
  CSTrainAppOptions opts_;
  // Count profiles for training
  profile_vector data_;
  // Count profiles for training
  auto_ptr< HMM<Alphabet> > hmm_;
  // Substitution matrix for pseudocount generation
  auto_ptr< SubstitutionMatrix<Alphabet> > subst_matrix_;
};  // CSTrainApp



template<class Alphabet>
void CSTrainApp<Alphabet>::parse_options(GetOpt_pp* options) {
  *options >> Option('i', "infile", opts_.infile, opts_.infile);
  *options >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  *options >> Option('d', "directory", opts_.directory, opts_.directory);
  *options >> Option('I', "informat", opts_.format, opts_.format);
  *options >> Option('M', "matchcol", opts_.matchcol_assignment,
                     opts_.matchcol_assignment);
  *options >> Option('K', "num-states", opts_.num_states, opts_.num_states);
  *options >> Option('W', "window-length", opts_.window_length,
                     opts_.window_length);
  *options >> Option('l', "likelihod-change", opts_.log_likelihood_change,
                     opts_.log_likelihood_change);
  *options >> Option('c', "connectivity", opts_.max_connectivity,
                     opts_.max_connectivity);
  *options >> Option('t', "transition-pc", opts_.transition_pseudocounts,
                     opts_.transition_pseudocounts);
  *options >> Option('s', "sample-rate", opts_.sample_rate,
                     opts_.sample_rate);
  *options >> Option('j', "jumpstart", opts_.hmmfile, opts_.hmmfile);
  *options >> Option('B', "blocks", opts_.num_blocks, opts_.num_blocks);
  *options >> Option('m', "matrix", opts_.blosum_type, opts_.blosum_type);
  *options >> Option('q', "mismatch-score", opts_.nucleotide_mismatch,
                     opts_.nucleotide_mismatch);
  *options >> Option('r', "match-score", opts_.nucleotide_match,
                     opts_.nucleotide_match);
  *options >> Option(' ', "data-pc", opts_.data_pseudocounts,
                     opts_.data_pseudocounts);
  *options >> Option(' ', "state-pc", opts_.state_pseudocounts,
                     opts_.state_pseudocounts);
  *options >> Option(' ', "min-scans", opts_.min_scans, opts_.min_scans);
  *options >> Option(' ', "max-scans", opts_.max_scans, opts_.max_scans);
  *options >> Option(' ', "weight-center", opts_.weight_center,
                     opts_.weight_center);
  *options >> Option(' ', "weight-decay", opts_.weight_decay,
                     opts_.weight_decay);
  *options >> Option(' ', "epsilon", opts_.epsilon_null, opts_.epsilon_null);
  *options >> Option(' ', "beta", opts_.beta, opts_.beta);
  *options >> OptionPresent(' ', "global-weights", opts_.global_weights);

  opts_.validate();

  if (!opts_.directory.empty() && *opts_.directory.rbegin() != kDirSep)
    opts_.directory.append(1, kDirSep);
  if (opts_.outfile.empty())
    opts_.outfile = opts_.directory +
      get_file_basename(opts_.infile, false) + "hmm";
  if (opts_.format == "auto")
    opts_.format = get_file_ext(opts_.infile);
}

template<class Alphabet>
void CSTrainApp<Alphabet>::print_description() const {
  fputs("Train a context HMM on dataset of full-length profiles, alignments, "
        "or sequences.\n", stream());
}

template<class Alphabet>
void CSTrainApp<Alphabet>::print_banner() const {
  fputs("Usage: cstrain -i <infile> -K <num_states> [options]\n", stream());
  fputs("       cstrain -i <infile> -j <hmmfile> [options]\n", stream());
}

template<class Alphabet>
void CSTrainApp<Alphabet>::print_substitution_matrix_options() const {
  fprintf(stream(), "  %-30s %s (def=%.0f)\n", "-q, --mismatch-score <int>",
          "Penalty for a nucleotide mismatch", opts_.nucleotide_mismatch);
  fprintf(stream(), "  %-30s %s (def=%.0f)\n",  "-r, --match-score <int>",
          "Reward for a nucleotide match", opts_.nucleotide_match);
}

template<>
void CSTrainApp<AminoAcid>::print_substitution_matrix_options() const {
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-m, --matrix <string>",
          "Substitution matrix: BLOSUM45, BLOSUM62, or BLOSUM80",
          opts_.blosum_type.c_str());
}

template<class Alphabet>
void CSTrainApp<Alphabet>::print_options() const {
  fprintf(stream(), "  %-30s %s\n", "-i, --infile <filename>",
          "Path to input file with training alignments or profiles");
  fprintf(stream(), "  %-30s %s\n", "-o, --outfile <filename>",
          "Path for output file with trained HMM");
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-d, --directory <directory>",
          "Directory for temporary and output files",
          opts_.directory.empty() ? "." : opts_.directory.c_str());
  fprintf(stream(), "  %-30s %s (def=%s)\n", "-I, --informat <string>",
          "Format of training data: prf, seq, fas, a2m, or a3m",
          opts_.format.c_str());
  fprintf(stream(), "  %-30s %s\n", "-M, --matchcol [0:100]",
          "Make all FASTA columns with less than X% gaps match columns");
  fprintf(stream(), "  %-30s %s\n", "", "(def: make columns with residue in "
          "first sequence match columns)");
  fprintf(stream(), "  %-30s %s\n", "-K, --num-states [0,inf[",
          "Number of states in the HMM to be trained");
  fprintf(stream(), "  %-30s %s (def=%i)\n", "-W, --window-length [0,inf[",
          "Length of context-window", opts_.window_length);
  fprintf(stream(), "  %-30s %s (def=%3.1g)\n", "-l, --likelihood [0,inf[",
          "Maximal likelihood change per column for convergence",
          opts_.log_likelihood_change);
  fprintf(stream(), "  %-30s %s (def=off)\n", "-c, --connectivity [1,K]",
          "Maximal state connectivity");
  fprintf(stream(), "  %-30s %s (def=%3.1f)\n", "-t, --transition-pc <float>",
          "Transition pseudocounts", opts_.transition_pseudocounts);
  fprintf(stream(), "  %-30s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]",
          "Fraction of profile windows sampled per subject",
          opts_.sample_rate);
  fprintf(stream(), "  %-30s %s\n", "-j, --jumpstart <filename>",
          "Jumpstart the HMM training with a serialized HMM.");
  fprintf(stream(), "  %-30s %s\n", "-B, --blocks [0,N]",
          "Number of blocks for online training (def: B=N^3/8)");

  print_substitution_matrix_options();

  fprintf(stream(), "  %-30s %s (def=%i)\n", "    --min-scans [0,inf[",
          "Minimal number of training data scans", opts_.min_scans);
  fprintf(stream(), "  %-30s %s (def=%i)\n", "    --max-scans [0,inf[",
          "Maximal number of training data scans", opts_.max_scans);
  fprintf(stream(), "  %-30s %s (def=%3.1f)\n", "    --state-pc [0,1]",
          "Pseudocounts for state profiles", opts_.state_pseudocounts);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --data-pc [0,1]",
          "Pseudocounts for training data", opts_.data_pseudocounts);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --weight-center [0,1]",
          "Weight of central profile column in context window",
          opts_.weight_center);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --weight-decay [0,1]",
          "Exponential decay of positional window weights",
          opts_.weight_decay);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --epsilon [0,1]",
          "Start value for learning rate epsilon in online training",
          opts_.epsilon_null);
  fprintf(stream(), "  %-30s %s (def=%4.2f)\n", "    --beta [0,1]",
          "Exponential decay of epsilon in online training",
          opts_.beta);
  fprintf(stream(), "  %-30s %s\n", "    --global-weights",
          "Use global instead of position-specific weights for profiles");
}

template<class Alphabet>
void CSTrainApp<Alphabet>::init_substitution_matrix() {
  subst_matrix_ =  auto_ptr< SubstitutionMatrix<Alphabet> >(
      new NucleotideMatrix(opts_.nucleotide_match,
                           opts_.nucleotide_mismatch));
}

template<>
void CSTrainApp<AminoAcid>::init_substitution_matrix() {
  BlosumMatrix::Type type = blosum_matrix_type_from_string(opts_.blosum_type);
  subst_matrix_ = auto_ptr< SubstitutionMatrix<AminoAcid> >(new BlosumMatrix(type));
}

template<class Alphabet>
void CSTrainApp<Alphabet>::read_training_data() {
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (!fin)
    throw Exception("Unable to read from input file '%s'!",
                    opts_.infile.c_str());

  if (opts_.format == "prf") {
    // Read data counts directly from serialized count profiles
    fprintf(stream(), "Reading training profiles from %s ...",
            get_file_basename(opts_.infile).c_str());
    fflush(stream());

    CountProfile<Alphabet>::readall(fin, &data_);
    fprintf(stream(), " %i profiles read\n", static_cast<int>(data_.size()));

  } else if (opts_.format == "seq") {
    // Read sequences and convert to counts
    fprintf(stream(), "Processing training sequences in %s ...\n",
            get_file_basename(opts_.infile).c_str());
    fflush(stream());

    int i = 0;
    while (!feof(fin)) {
      Sequence<Alphabet> seq(fin);
      shared_ptr< CountProfile<Alphabet> > cp_ptr(
          new CountProfile<Alphabet>(seq));
      data_.push_back(cp_ptr);

      i += 1;
      if (i % 2 == 0) {
        fputc('.', stream());
        fflush(stream());
      }
      if (i % 100 == 0) fprintf(stream(), " %i\n", i);
    }
    if (i % 100 != 0)
      fprintf(stream(), "%s %i\n",
              string(50 - iround((i % 100) / 2), ' ').c_str(), i);

  } else {
    // Read alignments and convert to counts
    fprintf(stream(), "Processing training alignments in %s ...\n",
            get_file_basename(opts_.infile).c_str());

    typename Alignment<Alphabet>::Format f =
      alignment_format_from_string<Alphabet>(opts_.format);
    int i = 0;
    while (!feof(fin)) {
      Alignment<Alphabet> ali(fin, f);
      if (f == Alignment<Alphabet>::FASTA) {
        if (opts_.matchcol_assignment < 0)
          ali.assign_match_columns_by_sequence();
        else
          ali.assign_match_columns_by_gap_rule(opts_.matchcol_assignment);
      }
      shared_ptr< CountProfile<Alphabet> > cp_ptr(
          new CountProfile<Alphabet>(ali, !opts_.global_weights));
      data_.push_back(cp_ptr);

      i += 1;
      if (i % 2 == 0) {
        fputc('.', stream());
        fflush(stream());
      }
      if (i % 100 == 0) fprintf(stream(), " %i\n", i);

      int c = fgetc(fin);
      if (c == EOF) break;
      ungetc(c, fin);
    }
    if (i % 100 != 0)
      fprintf(stream(), "%s %i\n",
              string(50 - iround((i % 100) / 2), ' ').c_str(), i);
  }

  fclose(fin);
}

template<class Alphabet>
void CSTrainApp<Alphabet>::init_hmm() {
  if (opts_.hmmfile.empty()) {  // seed the HMM from training data
    fprintf(stream(), "Initializing HMM by sampling %i context profiles from "
            "training profiles ...", opts_.num_states);
    fflush(stream());

    MatrixPseudocounts<Alphabet> pc(subst_matrix_.get());
    SamplingStateInitializer<Alphabet> state_init(data_,
                                                  opts_.sample_rate,
                                                  &pc,
                                                  opts_.state_pseudocounts);
    HomogeneousTransitionInitializer<Alphabet> transition_init;
    hmm_ = auto_ptr< HMM<Alphabet> >(new HMM<Alphabet>(opts_.num_states,
                                                       opts_.window_length,
                                                       state_init,
                                                       transition_init));
    hmm_->transform_states_to_logspace();
    fputc('\n', stream());

  } else {  // read HMM from jumpstart file
    FILE* fin = fopen(opts_.hmmfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read from hmmfile '%s'!",
                      opts_.hmmfile.c_str());
    fprintf(stream(), "Reading HMM from %s ...",
            get_file_basename(opts_.hmmfile).c_str());
    fflush(stream());

    hmm_ = auto_ptr< HMM<Alphabet> >(new HMM<Alphabet>(fin));

    fputc('\n', stream());
    fclose(fin);
  }
}

template<class Alphabet>
int CSTrainApp<Alphabet>::run() {
  init_substitution_matrix();
  read_training_data();
  init_hmm();

  // Add pseudocounts to training data
  fprintf(stream(), "Adding pseudocounts to training profiles (admix=%-.2f) ...",
          opts_.data_pseudocounts);
  fflush(stream());
  int num_data_cols = 0;
  MatrixPseudocounts<Alphabet> pc(subst_matrix_.get());
  for (profile_iterator ci = data_.begin(); ci != data_.end(); ++ci) {
    pc.add_to_profile(ConstantAdmixture(opts_.data_pseudocounts), ci->get());
    num_data_cols += (*ci)->num_cols();
  }
  fputc('\n', stream());

  // Run Baum-Welch training on HMM
  fprintf(stream(), "Running Baum-Welch training (K=%i, W=%i, N=%i, C=%i) ...",
          hmm_->num_states(), hmm_->num_cols(), static_cast<int>(data_.size()),
          num_data_cols);
  fflush(stream());
  fputs("\n\n", stream());
  BaumWelchTraining<Alphabet, CountProfile> bw(opts_, data_, *hmm_, stream());
  bw.run();

  // Write HMM to outfile
  FILE* fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout)
    throw Exception("Unable to write profiles to output file '%s'!",
                    opts_.outfile.c_str());
  hmm_->write(fout);
  fclose(fout);
  fprintf(stream(), "\nWrote HMM to %s\n", opts_.outfile.c_str());

  return 0;
}

}  // namespace cs
