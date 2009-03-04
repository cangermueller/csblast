#ifndef CS_CSTRAIN_H
#define CS_CSTRAIN_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation for HMM training.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "alignment.h"
#include "baum_welch_training.h"
#include "blosum_matrix.h"
#include "nucleotide_matrix.h"
#include "counts_profile.h"
#include "hmm.h"
#include "exception.h"
#include "forward_backward_algorithm.h"
#include "getopt_pp.h"
#include "matrix_pseudocounts.h"
#include "shared_ptr.h"
#include "util.h"

using namespace GetOpt;

namespace cs
{

template<class Alphabet_T>
class CSTrain : public BaumWelchParams,
                public SamplingStateInitializerParams
{
  public:
    // Constructs an algorithm object with parameters set by command line args.
    CSTrain();
    virtual ~CSTrain() { }

    // Runs the HMM training with the provided parameter settings.
    void run(GetOpt_pp& options, std::ostream& out = std::cout);
    // Prints usage message.
    std::ostream& usage(std::ostream& out);

  private:
    typedef std::vector< shared_ptr< CountsProfile<Alphabet_T> > > counts_vector;
    typedef typename counts_vector::iterator counts_iterator;

#ifdef _WIN32
    static const char DIR_SEP = '\\';
#else
    static const char DIR_SEP = '/';
#endif

    // Checks if all parameters are valid for running the training.
    void check();
    // Parses algorithm parameters from the command line options.
    void parse(GetOpt_pp& options);
    // Fills the vector with training data counts from input alignments or profiles.
    void read_training_data(counts_vector& v, std::ostream& out);
    // Returns a shared pointer to a substitution matrix.
    shared_ptr< SubstitutionMatrix<Alphabet_T> > get_substitution_matrix();
    // Prints usage for substitution matrix options.
    std::ostream& substitution_matrix_options(std::ostream& out);

    // The input alignment file with training data.
    std::string infile_;
    // The output file for the trained HMM.
    std::string outfile_;
    // Directory for output and temporary files
    std::string directory_;
    // File format of input alignment
    std::string format_;
    // HMM input file for restarting
    std::string hmmfile_;
    // Match column assignment for FASTA alignments
    int matchcol_assignment_;
    // The number of states in the HMM to train.
    int num_states_;
    // The number of columns in the context window.
    int window_length_;
    // Pseudocounts to be added to observed data counts.
    float data_pseudocounts_;
    // Use global instead of position specific weights for profile construction.
    bool global_weights_;
    // BLOSUM matrix for pseudocount generation.
    std::string blosum_type_;
    // Reward for a nucleotide match.
    float nucleotide_match_;
    // Penalty for a nucleotide mismatch
    float nucleotide_mismatch_;
    // The reporting level for logging.
    int log_level_;
};



template<class Alphabet_T>
CSTrain<Alphabet_T>::CSTrain()
        : format_("auto"),
          matchcol_assignment_(-1),
          num_states_(0),
          window_length_(1),
          data_pseudocounts_(0.01f),
          global_weights_(false),
          blosum_type_("BLOSUM62"),
          nucleotide_match_(1.0f),
          nucleotide_mismatch_(-3.0f),
          log_level_(Log::reporting_level())
{ }

template<class Alphabet_T>
void CSTrain<Alphabet_T>::run(GetOpt_pp& options, std::ostream& out)
{
    parse(options);
    shared_ptr< HMM<Alphabet_T> > hmm_ptr;
    shared_ptr< SubstitutionMatrix<Alphabet_T> > sm_ptr(get_substitution_matrix());
    MatrixPseudocounts<Alphabet_T> sm_pc(sm_ptr.get());
    counts_vector data;

    read_training_data(data, out);

    // construct or read the HMM
    if (hmmfile_.empty()) {
        out << strprintf("Initializing HMM by sampling %i context profiles from training profiles ...", num_states_);
        out.flush();
        LOG(INFO) << strprintf("Initializing HMM by sampling %i context profiles from training profiles ...", num_states_);
        SamplingStateInitializer<Alphabet_T> state_init(data, &sm_pc, *this);
        HomogeneousTransitionInitializer<Alphabet_T> transition_init;
        hmm_ptr = shared_ptr< HMM<Alphabet_T> >(new HMM<Alphabet_T>(num_states_, window_length_, state_init, transition_init));
        hmm_ptr->transform_states_to_logspace();
        out << std::endl;
    } else {
        std::ifstream fin(hmmfile_.c_str());
        out << strprintf("Reading HMM from %s ...", get_file_basename(hmmfile_).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading HMM from %s ...", get_file_basename(hmmfile_).c_str());
        hmm_ptr = shared_ptr< HMM<Alphabet_T> >(new HMM<Alphabet_T>(fin));
        out << std::endl;
    }

    // add pseudocounts to training data
    out << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", data_pseudocounts_);
    out.flush();
    LOG(INFO) << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", data_pseudocounts_);
    int num_data_cols = 0;
    for (counts_iterator ci = data.begin(); ci != data.end(); ++ci) {
        sm_pc.add_to_profile(**ci, ConstantAdmixture(data_pseudocounts_));
        (*ci)->convert_to_counts();
        num_data_cols += (*ci)->num_cols();
    }
    out << std::endl;

    // run Baum-Welch training on HMM
    out << strprintf("Running Baum-Welch training on HMM (K=%i, W=%i, N=%i) ...",
                     hmm_ptr->num_states(), hmm_ptr->num_cols(), num_data_cols);
    out.flush();
    LOG(INFO) << strprintf("Running Baum-Welch training on HMM (K=%i, W=%i, N=%i) ...",
                           hmm_ptr->num_states(), hmm_ptr->num_cols(), num_data_cols);
    out << std::endl << std::endl;
    TrainingProgressInfo<Alphabet_T> prg_info(*hmm_ptr, out);
    BaumWelchTraining<Alphabet_T, CountsProfile> training(*this);
    training.run(*hmm_ptr, data, &prg_info);

    // write HMM to outfile
    std::fstream fout(outfile_.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write HMM to output file '%s'!", outfile_.c_str());
    hmm_ptr->write(fout);
    fout.close();
    out << std::endl << strprintf("Wrote HMM to %s", outfile_.c_str()) << std::endl;
    LOG(INFO) << strprintf("Wrote HMM to %s", outfile_.c_str());
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::read_training_data(counts_vector& v, std::ostream& out)
{
    std::ifstream fin(infile_.c_str());

    if (format_ == "prf") {
        // read data counts directly from serialized counts profiles
        out << strprintf("Reading training profiles from %s ...", get_file_basename(infile_).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading training profiles from %s ...", get_file_basename(infile_).c_str());
        CountsProfile<Alphabet_T>::readall(fin, v);
        out << strprintf(" %i profiles read", v.size()) << std::endl;
        LOG(INFO) << strprintf("%i profiles read", v.size());

    } else if (format_ == "seq") {
        // read sequences and convert to counts
        out << strprintf("Processing training sequences in %s ...\n", get_file_basename(infile_).c_str());
        out.flush();
        LOG(INFO) << strprintf("Processing training sequences in %s ...", get_file_basename(infile_).c_str());

        int i = 0;
        while (fin.peek() && fin.good()) {
            Sequence<Alphabet_T> seq(fin);
            shared_ptr< CountsProfile<Alphabet_T> > cp_ptr(new CountsProfile<Alphabet_T>(seq));
            v.push_back(cp_ptr);

            i += 1;
            if (i % 2 == 0) {
                out << '.';
                out.flush();
            }
            if (i % 100 == 0) out << " " << i << std::endl;
        }
        if (i % 100 != 0) out << std::string(50 - iround((i % 100) / 2), ' ') << " " << i << std::endl;

    } else {
        // read alignments and convert to counts
        out << strprintf("Processing training alignments in %s ...\n", get_file_basename(infile_).c_str());
        LOG(INFO) << strprintf("Processing training alignments in %s ...", get_file_basename(infile_).c_str());

        typename Alignment<Alphabet_T>::Format f = alignment_format_from_string<Alphabet_T>(get_file_ext(infile_));
        int i = 0;
        while (fin.peek() && fin.good()) {
            Alignment<Alphabet_T> ali(fin, f);
            if (f == Alignment<Alphabet_T>::FASTA) {
                if (matchcol_assignment_ < 0)
                    ali.assign_match_columns_by_sequence();
                else
                    ali.assign_match_columns_by_gap_rule(matchcol_assignment_);
            }
            shared_ptr< CountsProfile<Alphabet_T> > cp_ptr(new CountsProfile<Alphabet_T>(ali, !global_weights_));
            v.push_back(cp_ptr);

            i += 1;
            if (i % 2 == 0) {
                out << '.';
                out.flush();
            }
            if (i % 100 == 0) out << " " << i << std::endl;
        }
        if (i % 100 != 0) out << std::string(50 - iround((i % 100) / 2), ' ') << " " << i << std::endl;
    }

    fin.close();
}

template<class Alphabet_T>
shared_ptr< SubstitutionMatrix<Alphabet_T> > CSTrain<Alphabet_T>::get_substitution_matrix()
{
    return shared_ptr< SubstitutionMatrix<Alphabet_T> >(new NucleotideMatrix(nucleotide_match_, nucleotide_mismatch_));
}

template<>
shared_ptr< SubstitutionMatrix<AminoAcid> > CSTrain<AminoAcid>::get_substitution_matrix()
{
    BlosumMatrix::Type type = blosum_matrix_type_from_string(blosum_type_);
    shared_ptr< SubstitutionMatrix<AminoAcid> > matrix_ptr(new BlosumMatrix(type));
    return matrix_ptr;
}

template<class Alphabet_T>
std::ostream& CSTrain<Alphabet_T>::usage(std::ostream& out)
{
    out << "Train an HMM of context-states on a dataset of alignments or sequences.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain -i <infile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-38s %s\n",             "-i, --infile <filename>", "Path to input file with training alignments or profiles");
    out << strprintf("  %-38s %s\n",             "-o, --outfile <filename>", "Path to output file with trained HMM");
    out << strprintf("  %-38s %s (def=%s)\n",    "-d, --directory <directory>", "Directory for temporary and output files",
                     directory_.empty() ? "." : directory_.c_str());
    out << strprintf("  %-38s %s (def=%s)\n",    "-f, --format <string>", "Format of training data: prf, seq, fas, a2m, or a3m",
                     format_.c_str());
    out << strprintf("  %-38s %s\n",             "-M, --matchcol-assignment [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    out << strprintf("  %-38s %s\n",             "", "(def: make columns with residue in first sequence match columns)");
    out << strprintf("  %-38s %s\n",             "-K, --num-states [0,inf[", "Number of states in the HMM to be trained");
    out << strprintf("  %-38s %s (def=%i)\n",    "-W, --window-length [0,inf[", "Length of context-window",
                     window_length_);
    out << strprintf("  %-38s %s (def=%3.1g)\n", "-l, --likelihood [0,inf[", "Maximal likelihood change per column for convergence",
                     log_likelihood_threshold_);
    out << strprintf("  %-38s %s (def=off)\n",   "-c, --connectivity [1,K]", "Maximal state connectivity",
                     max_connectivity_);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "-t, --transition-pseudocounts <float>", "Transition pseudocounts",
                     transition_pseudocounts_);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]", "Context window sample rate for initialization",
                     sample_rate_);
    out << strprintf("  %-38s %s\n",             "-j, --jumpstart <filename>", "Jumpstart the HMM training with a serialized HMM.");
    substitution_matrix_options(out);
    out << strprintf("  %-38s %s (def=%i)\n",    "    --min-iterations [0,inf[", "Minimal number of training iterations",
                     min_iterations_);
    out << strprintf("  %-38s %s (def=%i)\n",    "    --max-iterations [0,inf[", "Maximal number of training iterations",
                     max_iterations_);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "    --state-pseudocounts [0,1]", "Pseudocounts for state profiles",
                     state_pseudocounts_);
    out << strprintf("  %-38s %s (def=%4.2f)\n", "    --data-pseudocounts [0,1]", "Pseudocounts for training data",
                     data_pseudocounts_);
    out << strprintf("  %-38s %s\n",             "    --global-weights", "Use global instead of position-specific weights for profiles");
    out << strprintf("  %-38s %s (def=%i)\n",    "    --log-max-level <int>", "Maximal reporting level for logging",
                     log_level_);

    return out;
}

template<class Alphabet_T>
std::ostream& CSTrain<Alphabet_T>::substitution_matrix_options(std::ostream& out)
{
    out << strprintf("  %-38s %s (def=%.0f)\n",  "-q, --mismatch-score <int>", "Penalty for a nucleotide mismatch",
                     nucleotide_mismatch_);
    out << strprintf("  %-38s %s (def=%.0f)\n",  "-r, --match-score <int>", "Reward for a nucleotide match",
                     nucleotide_match_);
    return out;
}

template<>
std::ostream& CSTrain<AminoAcid>::substitution_matrix_options(std::ostream& out)
{
    out << strprintf("  %-38s %s (def=%s)\n",    "-m, --matrix <string>", "Substitution matrix: BLOSUM45, BLOSUM62, or BLOSUM80",
                     blosum_type_.c_str());
    return out;
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::parse(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile_, infile_);
    options >> Option('o', "outfile", outfile_, outfile_);
    options >> Option('d', "directory", directory_, directory_);
    options >> Option('f', "format", format_, format_);
    options >> Option('M', "matchcol-assignment", matchcol_assignment_, matchcol_assignment_);
    options >> Option('K', "num-states", num_states_, num_states_);
    options >> Option('W', "window-length", window_length_, window_length_);
    options >> Option('l', "likelihod-change", log_likelihood_threshold_, log_likelihood_threshold_);
    options >> Option('c', "connectivity", max_connectivity_, max_connectivity_);
    options >> Option('t', "transition_pseudocounts", transition_pseudocounts_, transition_pseudocounts_);
    options >> Option('s', "sample-rate", sample_rate_, sample_rate_);
    options >> Option('j', "jumpstart", hmmfile_, hmmfile_);
    options >> Option('m', "matrix", blosum_type_, blosum_type_);
    options >> Option('q', "mismatch-score", nucleotide_mismatch_, nucleotide_mismatch_);
    options >> Option('r', "match-score", nucleotide_match_, nucleotide_match_);
    options >> Option(' ', "data-pseudocounts", data_pseudocounts_, data_pseudocounts_);
    options >> Option(' ', "state-pseudocounts", state_pseudocounts_, state_pseudocounts_);
    options >> Option(' ', "min-iterations", min_iterations_, min_iterations_);
    options >> Option(' ', "max-iterations", max_iterations_, max_iterations_);
    options >> OptionPresent(' ', "global-weights", global_weights_);
    options >> Option(' ', "max-log-level", log_level_, log_level_);
    Log::reporting_level() = Log::from_integer(log_level_);

    check();
    if (!directory_.empty() && *directory_.rbegin() != DIR_SEP) directory_.append(1, DIR_SEP);
    if (outfile_.empty()) outfile_ = directory_ + get_file_basename(infile_, false) + "hmm";
    if (format_ == "auto") format_ = get_file_ext(infile_);
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::check()
{
    if (num_states_ == 0 && hmmfile_.empty()) throw cs::Exception("No value for number of HMM states provided!");
    if (infile_.empty()) throw cs::Exception("No input file with training data provided!");
}

}  // cs

#endif
