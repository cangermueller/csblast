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
    std::ifstream ali_in(infile_.c_str());
    ali_vector alis = Alignment<Alphabet_T>::readall(ali_in, Alignment<Alphabet_T>::FASTA);
    ali_in.close();
    for (ali_iterator ai = alis.begin(); ai != alis.end(); ++ai)
        (*ai)->assign_match_columns_by_sequence();

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

    // Add pseudocounts to observed counts
    for (counts_iterator ci = counts.begin(); ci != counts.end(); ++ci) {
        sm_pc.add_to_profile(**ci, ConstantAdmixture(data_pseudocounts_));
        (*ci)->convert_to_counts();
    }

    // Perform Baum-Welch training on HMM
    TrainingProgressInfo<Alphabet_T> prg_info(hmm);
    BaumWelchTraining<Alphabet_T, CountsProfile> training(*this);
    training.run(hmm, counts, &prg_info);
}

template<class Alphabet_T>
std::ostream& CSTrain<Alphabet_T>::usage(std::ostream& out)
{
    out << "Train an HMM of context-states on a dataset of alignments.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain_aa -i <alifile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << cs::strprintf("  %-30s %s\n",          "-K, --num-states <int>", "Number of states in the HMM to be trained");
    out << cs::strprintf("  %-30s %s (def=%i)\n", "-l, --log-max-level <int>", "Maximal reporting level for logging", log_level_);

    return out;
}

template<class Alphabet_T>
void CSTrain<Alphabet_T>::parse(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile_, infile_);
    options >> Option('o', "outfile", outfile_, outfile_);
    options >> Option('K', "num-states", num_states_, num_states_);
    options >> Option('W', "window-length", window_length_, window_length_);
    options >> Option('d', "data-pseudocounts", data_pseudocounts_, data_pseudocounts_);
    options >> Option('s', "state-pseudocounts", state_pseudocounts_, state_pseudocounts_);
    options >> Option('r', "sample-rate", sample_rate_, sample_rate_);
    options >> Option('a', "alpha", alpha_, alpha_);
    options >> Option('l', "max-log-level", log_level_, log_level_);
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
