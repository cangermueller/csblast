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

struct Params : public BaumWelchParams
{
    Params()
            : format("auto"),
              matchcol_assignment(-1),
              num_states(0),
              window_length(1),
              sample_rate(0.2f),
              state_pseudocounts(1.0f),
              data_pseudocounts(0.01f),
              global_weights(false),
              blosum_type("BLOSUM62"),
              nucleotide_match(1.0f),
              nucleotide_mismatch(-3.0f),
              log_level(Log::reporting_level())
    { }

    virtual ~Params() { }

    // Checks if all parameters are valid for running the training.
    void check();
    // Parses command line arguments and sets parameters accordingly.
    void parse_options(GetOpt_pp& options);

#ifdef _WIN32
    static const char DIR_SEP = '\\';
#else
    static const char DIR_SEP = '/';
#endif

    // The input alignment file with training data.
    std::string infile;
    // The output file for the trained HMM.
    std::string outfile;
    // Directory for output and temporary files
    std::string directory;
    // File format of input alignment
    std::string format;
    // HMM input file for restarting
    std::string hmmfile;
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
    std::string blosum_type;
    // Reward for a nucleotide match.
    float nucleotide_match;
    // Penalty for a nucleotide mismatch
    float nucleotide_mismatch;
    // The reporting level for logging.
    int log_level;
};



void Params::parse_options(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile, infile);
    options >> Option('o', "outfile", outfile, outfile);
    options >> Option('d', "directory", directory, directory);
    options >> Option('f', "format", format, format);
    options >> Option('M', "matchcol-assignment", matchcol_assignment, matchcol_assignment);
    options >> Option('K', "num-states", num_states, num_states);
    options >> Option('W', "window-length", window_length, window_length);
    options >> Option('l', "likelihod-change", log_likelihood_threshold, log_likelihood_threshold);
    options >> Option('c', "connectivity", max_connectivity, max_connectivity);
    options >> Option('t', "transition_pseudocounts", transition_pseudocounts, transition_pseudocounts);
    options >> Option('s', "sample-rate", sample_rate, sample_rate);
    options >> Option('j', "jumpstart", hmmfile, hmmfile);
    options >> Option('m', "matrix", blosum_type, blosum_type);
    options >> Option('q', "mismatch-score", nucleotide_mismatch, nucleotide_mismatch);
    options >> Option('r', "match-score", nucleotide_match, nucleotide_match);
    options >> Option(' ', "data-pseudocounts", data_pseudocounts, data_pseudocounts);
    options >> Option(' ', "state-pseudocounts", state_pseudocounts, state_pseudocounts);
    options >> Option(' ', "min-iterations", min_iterations, min_iterations);
    options >> Option(' ', "max-iterations", max_iterations, max_iterations);
    options >> Option(' ', "weight-center", weight_center, weight_center);
    options >> Option(' ', "weight-decay", weight_decay, weight_decay);
    options >> OptionPresent(' ', "global-weights", global_weights);
    options >> Option(' ', "log-level", log_level, log_level);
    Log::reporting_level() = Log::from_integer(log_level);

    check();
    if (!directory.empty() && *directory.rbegin() != DIR_SEP) directory.append(1, DIR_SEP);
    if (outfile.empty()) outfile = directory + get_file_basename(infile, false) + "hmm";
    if (format == "auto") format = get_file_ext(infile);
}

void Params::check()
{
    if (num_states == 0 && hmmfile.empty()) throw Exception("No value for number of HMM states provided!");
    if (infile.empty()) throw Exception("No input file with training data provided!");
}



template<class Alphabet_T>
std::ostream& usage(const Params& params, std::ostream& out = std::cout);

template<class Alphabet_T>
std::ostream& substitution_matrix_options(const Params& params, std::ostream& out);

template<>
std::ostream& substitution_matrix_options<AminoAcid>(const Params& params, std::ostream& out);

template<class Alphabet_T>
std::ostream& usage(const Params& params, std::ostream& out = std::cout)
{
    out << "Train an HMM of context-states on a dataset of alignments or sequences.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain -i <infile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-38s %s\n",             "-i, --infile <filename>", "Path to input file with training alignments or profiles");
    out << strprintf("  %-38s %s\n",             "-o, --outfile <filename>", "Path for output file with trained HMM");
    out << strprintf("  %-38s %s (def=%s)\n",    "-d, --directory <directory>", "Directory for temporary and output files",
                     params.directory.empty() ? "." : params.directory.c_str());
    out << strprintf("  %-38s %s (def=%s)\n",    "-f, --format <string>", "Format of training data: prf, seq, fas, a2m, or a3m",
                     params.format.c_str());
    out << strprintf("  %-38s %s\n",             "-M, --matchcol-assignment [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    out << strprintf("  %-38s %s\n",             "", "(def: make columns with residue in first sequence match columns)");
    out << strprintf("  %-38s %s\n",             "-K, --num-states [0,inf[", "Number of states in the HMM to be trained");
    out << strprintf("  %-38s %s (def=%i)\n",    "-W, --window-length [0,inf[", "Length of context-window",
                     params.window_length);
    out << strprintf("  %-38s %s (def=%3.1g)\n", "-l, --likelihood [0,inf[", "Maximal likelihood change per column for convergence",
                     params.log_likelihood_threshold);
    out << strprintf("  %-38s %s (def=off)\n",   "-c, --connectivity [1,K]", "Maximal state connectivity",
                     params.max_connectivity);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "-t, --transition-pseudocounts <float>", "Transition pseudocounts",
                     params.transition_pseudocounts);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]", "Fraction of profile windows sampled per subject for HMM initialization",
                     params.sample_rate);
    out << strprintf("  %-38s %s\n",             "-j, --jumpstart <filename>", "Jumpstart the HMM training with a serialized HMM.");

    substitution_matrix_options<Alphabet_T>(params, out);

    out << strprintf("  %-38s %s (def=%i)\n",    "    --min-iterations [0,inf[", "Minimal number of training iterations",
                     params.min_iterations);
    out << strprintf("  %-38s %s (def=%i)\n",    "    --max-iterations [0,inf[", "Maximal number of training iterations",
                     params.max_iterations);
    out << strprintf("  %-38s %s (def=%3.1f)\n", "    --state-pseudocounts [0,1]", "Pseudocounts for state profiles",
                     params.state_pseudocounts);
    out << strprintf("  %-38s %s (def=%4.2f)\n", "    --data-pseudocounts [0,1]", "Pseudocounts for training data",
                     params.data_pseudocounts);
    out << strprintf("  %-38s %s (def=%4.2f)\n", "    --weight-center [0,1]", "Weight of central profile column in context window.",
                     params.weight_center);
    out << strprintf("  %-38s %s (def=%4.2f)\n", "    --weight-decay [0,1]", "Exponential decay of positional window weights.",
                     params.weight_decay);
    out << strprintf("  %-38s %s\n",             "    --global-weights", "Use global instead of position-specific weights for profiles");
    out << strprintf("  %-38s %s (def=%i)\n",    "    --log-level <int>", "Maximal reporting level for logging",
                     params.log_level);

    return out;
}

template<class Alphabet_T>
std::ostream& substitution_matrix_options(const Params& params, std::ostream& out)
{
    out << strprintf("  %-38s %s (def=%.0f)\n",  "-q, --mismatch-score <int>", "Penalty for a nucleotide mismatch",
                     params.nucleotide_mismatch);
    out << strprintf("  %-38s %s (def=%.0f)\n",  "-r, --match-score <int>", "Reward for a nucleotide match",
                     params.nucleotide_match);
    return out;
}

template<>
std::ostream& substitution_matrix_options<AminoAcid>(const Params& params, std::ostream& out)
{
    out << strprintf("  %-38s %s (def=%s)\n",    "-m, --matrix <string>", "Substitution matrix: BLOSUM45, BLOSUM62, or BLOSUM80",
                     params.blosum_type.c_str());
    return out;
}



template<class Alphabet_T>
void cstrain(const Params& params, std::ostream& out = std::cout);

template<class Alphabet_T>
void read_training_data(const Params& params, std::vector< shared_ptr< CountsProfile<Alphabet_T> > >& v, std::ostream& out);

template<class Alphabet_T>
shared_ptr< SubstitutionMatrix<Alphabet_T> > get_substitution_matrix(const Params& params);

template<>
shared_ptr< SubstitutionMatrix<AminoAcid> > get_substitution_matrix<AminoAcid>(const Params& params);

template<class Alphabet_T>
void cstrain(const Params& params, std::ostream& out)
{
    typedef std::vector< shared_ptr< CountsProfile<Alphabet_T> > > counts_vector;
    typedef typename counts_vector::iterator counts_iterator;

    shared_ptr< HMM<Alphabet_T> > hmm_ptr;
    shared_ptr< SubstitutionMatrix<Alphabet_T> > sm_ptr(get_substitution_matrix<Alphabet_T>(params));
    MatrixPseudocounts<Alphabet_T> sm_pc(sm_ptr.get());
    counts_vector data;

    read_training_data(params, data, out);

    // construct or read the HMM
    if (params.hmmfile.empty()) {
        out << strprintf("Initializing HMM by sampling %i context profiles from training profiles ...", params.num_states);
        out.flush();
        LOG(INFO) << strprintf("Initializing HMM by sampling %i context profiles from training profiles ...", params.num_states);
        SamplingStateInitializer<Alphabet_T> state_init(data, params.sample_rate, &sm_pc, params.state_pseudocounts);
        HomogeneousTransitionInitializer<Alphabet_T> transition_init;
        hmm_ptr = shared_ptr< HMM<Alphabet_T> >(new HMM<Alphabet_T>(params.num_states, params.window_length, state_init, transition_init));
        hmm_ptr->transform_states_to_logspace();
        out << std::endl;
    } else {
        std::ifstream fin(params.hmmfile.c_str());
        out << strprintf("Reading HMM from %s ...", get_file_basename(params.hmmfile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading HMM from %s ...", get_file_basename(params.hmmfile).c_str());
        hmm_ptr = shared_ptr< HMM<Alphabet_T> >(new HMM<Alphabet_T>(fin));
        out << std::endl;
    }

    // add pseudocounts to training data
    out << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", params.data_pseudocounts);
    out.flush();
    LOG(INFO) << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", params.data_pseudocounts);
    int num_data_cols = 0;
    for (counts_iterator ci = data.begin(); ci != data.end(); ++ci) {
        sm_pc.add_to_profile(**ci, ConstantAdmixture(params.data_pseudocounts));
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
    BaumWelchTraining<Alphabet_T, CountsProfile> training(*hmm_ptr, data, params, &prg_info);
    training.run();

    // write HMM to outfile
    std::ofstream fout(params.outfile.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write HMM to output file '%s'!", params.outfile.c_str());
    hmm_ptr->write(fout);
    fout.close();
    out << std::endl << strprintf("Wrote HMM to %s", params.outfile.c_str()) << std::endl;
    LOG(INFO) << strprintf("Wrote HMM to %s", params.outfile.c_str());
}

template<class Alphabet_T>
void read_training_data(const Params& params, std::vector< shared_ptr< CountsProfile<Alphabet_T> > >& v, std::ostream& out)
{
    std::ifstream fin(params.infile.c_str());
    if (!fin) throw Exception("Unable to read from input file '%s'!", params.infile.c_str());

    if (params.format == "prf") {
        // read data counts directly from serialized counts profiles
        out << strprintf("Reading training profiles from %s ...", get_file_basename(params.infile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading training profiles from %s ...", get_file_basename(params.infile).c_str());
        CountsProfile<Alphabet_T>::readall(fin, v);
        out << strprintf(" %i profiles read", v.size()) << std::endl;
        LOG(INFO) << strprintf("%i profiles read", v.size());

    } else if (params.format == "seq") {
        // read sequences and convert to counts
        out << strprintf("Processing training sequences in %s ...\n", get_file_basename(params.infile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Processing training sequences in %s ...", get_file_basename(params.infile).c_str());

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
        out << strprintf("Processing training alignments in %s ...\n", get_file_basename(params.infile).c_str());
        LOG(INFO) << strprintf("Processing training alignments in %s ...", get_file_basename(params.infile).c_str());

        typename Alignment<Alphabet_T>::Format f = alignment_format_from_string<Alphabet_T>(params.format);
        int i = 0;
        while (fin.peek() && fin.good()) {
            Alignment<Alphabet_T> ali(fin, f);
            if (f == Alignment<Alphabet_T>::FASTA) {
                if (params.matchcol_assignment < 0)
                    ali.assign_match_columns_by_sequence();
                else
                    ali.assign_match_columns_by_gap_rule(params.matchcol_assignment);
            }
            shared_ptr< CountsProfile<Alphabet_T> > cp_ptr(new CountsProfile<Alphabet_T>(ali, !params.global_weights));
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
shared_ptr< SubstitutionMatrix<Alphabet_T> > get_substitution_matrix(const Params& params)
{
    return shared_ptr< SubstitutionMatrix<Alphabet_T> >(new NucleotideMatrix(params.nucleotide_match, params.nucleotide_mismatch));
}

template<>
shared_ptr< SubstitutionMatrix<AminoAcid> > get_substitution_matrix<AminoAcid>(const Params& params)
{
    BlosumMatrix::Type type = blosum_matrix_type_from_string(params.blosum_type);
    shared_ptr< SubstitutionMatrix<AminoAcid> > matrix_ptr(new BlosumMatrix(type));
    return matrix_ptr;
}


}  // cs

#endif
