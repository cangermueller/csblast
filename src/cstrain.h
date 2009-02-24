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
    void run(GetOpt_pp& options);
    // Prints usage message.
    std::ostream& usage(std::ostream& out);

  private:
    // Checks if all parameters are valid for running the training.
    void check();
    // Parses algorithm parameters from the command line options.
    void parse(GetOpt_pp& options);

    // The input alignment file with training data.
    std::string infile_;
    // The output file for the trained HMM.
    std::string outfile_;
    // The number of states in the HMM to train.
    int num_states_;
    // The number of columns in the context window.
    int window_length_;
    // Pseudocounts to be added to observed data counts.
    float data_pseudocounts_;
    // The reporting level for logging.
    int log_level_;
};



template<class Alphabet_T>
CSTrain<Alphabet_T>::CSTrain()
        : num_states_(0),
          window_length_(1),
          data_pseudocounts_(0.01f),
          log_level_(Log::reporting_level())
{ }

template<class Alphabet_T>
void CSTrain<Alphabet_T>::run(GetOpt_pp& options)
{
    parse(options);
    typedef std::vector< shared_ptr< Alignment<Alphabet_T> > > ali_vector;
    typedef typename std::vector< shared_ptr< Alignment<Alphabet_T> > >::iterator ali_iterator;
    typedef std::vector< shared_ptr< CountsProfile<Alphabet_T> > > counts_vector;
    typedef typename std::vector< shared_ptr<CountsProfile<Alphabet_T> > >::iterator counts_iterator;

    // Read input alignments
    std::ifstream fin(infile_.c_str());
    ali_vector alis = Alignment<Alphabet_T>::readall(fin, Alignment<Alphabet_T>::A3M);
    fin.close();
    // for (ali_iterator ai = alis.begin(); ai != alis.end(); ++ai)
    //     (*ai)->assign_match_columns_by_sequence();

    // Convert alignments to counts profiles
    counts_vector counts;
    for (ali_iterator ai = alis.begin(); ai != alis.end(); ++ai) {
        shared_ptr< CountsProfile<Alphabet_T> > cp_ptr(new CountsProfile<Alphabet_T>(**ai, true));
        counts.push_back(cp_ptr);
    }
    alis.clear();

    // Build the HMM
    BlosumMatrix sm;
    MatrixPseudocounts<AminoAcid> sm_pc(&sm);
    SamplingStateInitializer<Alphabet_T> state_init(counts, &sm_pc, *this);
    HomogeneousTransitionInitializer<Alphabet_T> transition_init;
    HMM<Alphabet_T> hmm(num_states_, window_length_, state_init, transition_init);
    hmm.transform_states_to_logspace();

    // Add pseudocounts to observed counts
    for (counts_iterator ci = counts.begin(); ci != counts.end(); ++ci) {
        sm_pc.add_to_profile(**ci, ConstantAdmixture(data_pseudocounts_));
        (*ci)->convert_to_counts();
    }

    // Perform Baum-Welch training on HMM
    TrainingProgressInfo<Alphabet_T> prg_info(hmm);
    BaumWelchTraining<Alphabet_T, CountsProfile> training(*this);
    training.run(hmm, counts, &prg_info);

    // Write HMM to outfile.
    std::fstream fout(outfile_.c_str(), std::ios_base::out);
    if(!fout)
        throw Exception("Unable to write to output file '%s'!", outfile_.c_str());
    else
        hmm.write(fout);
    fout.close();
}

template<class Alphabet_T>
std::ostream& CSTrain<Alphabet_T>::usage(std::ostream& out)
{
    out << "Train an HMM of context-states on a dataset of alignments.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain -i <infile> -o <outfile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-30s %s\n",             "-i, --infile <filename>", "Input file with alignments for training");
    out << strprintf("  %-30s %s\n",             "-o, --outfile <filename>", "Output file with trained HMM");
    out << strprintf("  %-30s %s\n",             "-K, --num-states [0,inf[]", "Number of states in the HMM to be trained");
    out << strprintf("  %-30s %s (def=%i)\n",    "-W, --window-length [0,inf[", "Length of context-window", window_length_);
    out << strprintf("  %-30s %s (def=%3.1f)\n", "-t, --threshold [0,inf[", "Likelihood change threshold for convergence", log_likelihood_threshold_);
    out << strprintf("  %-30s %s (def=%3.1f)\n", "-a, --alpha <float>", "Parameter alpha for transitions prior", alpha_);
    out << strprintf("  %-30s %s (def=%3.1f)\n", "-r, --sample-rate [0,1]", "Context window sample rate in HMM initialization", sample_rate_);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --min-iterations [0,inf[", "Minimal number of training iterations", min_iterations_);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --max-iterations [0,inf[", "Maximal number of training iterations", max_iterations_);
    out << strprintf("  %-30s %s (def=%3.1f)\n", "    --state-pseudocounts [0,1]", "Pseudocounts for state profiles", state_pseudocounts_);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --data-pseudocounts [0,1]", "Pseudocounts for data counts", data_pseudocounts_);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --log-max-level <int>", "Maximal reporting level for logging", log_level_);

    return out;
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::parse(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile_, infile_);
    options >> Option('o', "outfile", outfile_, outfile_);
    options >> Option('K', "num-states", num_states_, num_states_);
    options >> Option('W', "window-length", window_length_, window_length_);
    options >> Option('t', "threshold", log_likelihood_threshold_, log_likelihood_threshold_);
    options >> Option('a', "alpha", alpha_, alpha_);
    options >> Option('r', "sample-rate", sample_rate_, sample_rate_);
    options >> Option(' ', "data-pseudocounts", data_pseudocounts_, data_pseudocounts_);
    options >> Option(' ', "state-pseudocounts", state_pseudocounts_, state_pseudocounts_);
    options >> Option(' ', "min-iterations", min_iterations_, min_iterations_);
    options >> Option(' ', "max-iterations", max_iterations_, max_iterations_);
    options >> Option(' ', "max-log-level", log_level_, log_level_);
    Log::reporting_level() = Log::from_integer(log_level_);

    check();
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::check()
{
    if (num_states_ == 0) throw cs::Exception("No value for number of HMM states provided!");
    if (infile_.empty()) throw cs::Exception("No input file provided!");
    if (outfile_.empty()) throw cs::Exception("No output file provided!");
}

}  // cs

#endif
