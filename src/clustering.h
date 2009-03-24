#ifndef CS_BAUM_WELCH_TRAINING_H
#define CS_BAUM_WELCH_TRAINING_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of Baum-Welch training for HMMs.

#include <cmath>

#include <iostream>
#include <vector>

#include "forward_backward_algorithm.h"
#include "log.h"
#include "context_profile.h"
#include "profile_matcher.h"
#include "profile_library.h"
#include "progress_bar.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

// Forward declaration;
template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class Clustering;



class ClusteringProgressTable : public ProgressBar
{
  public:
    ClusteringProgressTable(std::ostream& out = std::cout)
            : ProgressBar(out, 30)
    { }

    // Prints header information.
    void print_header()
    {
        out_ << strprintf("%-4s %4s %4s %7s  %-30s  %9s  %8s\n",
                          "Scan", "Itrs", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
        out_ << std::string(82, '-') << std::endl;
    }

    // Starts a new row and prints statistics to outstream.
    void print_row_begin(int scan, int iter, int blocks, float epsilon)
    {
        reset();
        out_ << strprintf("%-4i %4i %4i %7.4f  ", scan, iter, blocks, epsilon);
        out_.flush();
    }

    // Ends the current row and prints likelihood.
    void print_row_end(float log_likelihood)
    {
        out_ << strprintf("  %9.5f\n", log_likelihood);
        out_.flush();
    }

    // Ends the current row and prints likelihood with change from previous scan.
    void print_row_end(float log_likelihood, float delta)
    {
        out_ << strprintf("  %9.5f  %+8.5f\n", log_likelihood, delta);
        out_.flush();
    }
};



struct ClusteringParams : public ForwardBackwardParams
{
    BaumWelchParams()
            : ForwardBackwardParams(),
              min_scans(50),
              max_scans(500),
              log_likelihood_change(2e-4f),
              num_blocks(0),
              epsilon_null(0.5f),
              beta(0.2f),
              epsilon_batch(0.05f)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_scans(params.min_scans),
              max_scans(params.max_scans),
              log_likelihood_change(params.log_likelihood_change),
              num_blocks(params.num_blocks),
              epsilon_null(params.epsilon_null),
              beta(params.beta),
              epsilon_batch(params.epsilon_batch)
    { }

    // Minimal number of data scans.
    int min_scans;
    // Maximum number of data scans.
    int max_scans;
    // Log-likelihood change per column for convergence.
    float log_likelihood_change;
    // Number of blocks into which the training data are divided (default:  B=N^(3/8)).
    int num_blocks;
    // Initial value for learning rate epsilon (1-epsilon is preserved fraction of sufficient statistics).
    float epsilon_null;
    // Paramter governing exponential decay of epsilon per scan.
    float beta;
    // Learning rate epsilon below which training is switched to batch mode.
    float epsilon_batch;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class Clustering
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
    typedef typename std::vector<data_vector> blocks_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;

    // Initializes a new training object.
    BaumWelchTraining(const BaumWelchParams& params,
                      const data_vector& data,
                      HMM<Alphabet_T>& hmm,
                      TrainingProgressInfo* progress_info = NULL);

    virtual ~BaumWelchTraining()
    { }

    // Trains the HMM with the data provided until one of the termination criterions is fullfilled.
    void run();

  private:
    // Runs forward backward algorithm on provided data.
    void expectation_step(const data_vector& block, float epsilon);
    // Calculates and assigns new HMM parameters by Maxmimum-Likelihood estimation.
    void maximization_step();
    // Adds the contribution of a subject's forward-backward matrices to prio probabilities of states.
    void add_contribution_to_priors(const ForwardBackwardMatrices& m);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Updates global sufficient statistics with sufficient statistics calculated on current block.
    void update_sufficient_statistics(float epsilon);
    // Fills the blocks vector with training data.
    void setup_blocks(bool batch_mode = false);
    // Prepares all members for HMM training.
    void init();
    // Calculates the change between the current log-likelihood and the previous-iteration log-likelihood.
    float log_likelihood_change();
    // Returns true if any termination conditions is fullfilled.
    bool terminate();

    // Parameter wrapper for Baum-Welch training.
    const BaumWelchParams& params_;
    // Training data (either sequences or counts profiles)
    const data_vector& data_;
    // HMM to be trained
    HMM<Alphabet_T>& hmm_;
    // Blocks of training data
    blocks_vector blocks_;
    // Progress info object for output of progress table
    TrainingProgressInfo* progress_info_;
    // Global expected sufficient statistics for transitions
    sparse_matrix<float> transition_stats_;
    // Global expeted sufficient statistics for emissions and state priors
    profiles_vector profile_stats_;
    // Expected sufficient statistics for transitions based on current block
    sparse_matrix<float> transition_stats_block_;
    // Expeted sufficient statistics for emissions and state priors based on current block
    profiles_vector profile_stats_block_;
    // Incremental likelihood in current iteration
    float log_likelihood_;
    // Complete likelihood after last complete scan of training data
    float log_likelihood_prev_;
    // Number of traning iterations performed so far.
    int iterations_;
    // Number of current data scan.
    int scan_;
    // Effective number of training columns.
    int num_eff_cols_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline BaumWelchTraining<Alphabet_T, Subject_T>::BaumWelchTraining(const BaumWelchParams& params,
                                                                   const data_vector& data,
                                                                   HMM<Alphabet_T>& hmm,
                                                                   TrainingProgressInfo* progress_info)
        : params_(params),
          data_(data),
          hmm_(hmm),
          blocks_(),
          progress_info_(progress_info),
          transition_stats_(hmm.num_states(), hmm.num_states()),
          profile_stats_(),
          transition_stats_block_(hmm.num_states(), hmm.num_states()),
          profile_stats_block_(),
          log_likelihood_(0.0f),
          log_likelihood_prev_(0.0f),
          iterations_(0),
          scan_(1),
          num_eff_cols_(0)
{
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::run()
{
    LOG(DEBUG) << "Running Baum-Welch training on ...";
    LOG(DEBUG) << hmm_;

    // Calculate log-likelihood baseline by batch EM without blocks
    float epsilon = 1.0f;
    if (progress_info_) progress_info_->start_scan(scan_, hmm_.iterations(), hmm_.connectivity(), 1, epsilon);
    expectation_step(data_, epsilon);
    maximization_step();
    ++iterations_;
    if (progress_info_) progress_info_->complete_scan(log_likelihood_);

    // Perform E-step and M-step on each training block until convergence
    setup_blocks();
    while (!terminate()) {
        if (static_cast<int>(blocks_.size()) > 1) epsilon = params_.epsilon_null * exp(-params_.beta * (scan_ - 1));
        log_likelihood_prev_ = log_likelihood_;
        log_likelihood_      = 0.0;
        ++scan_;

        if (progress_info_) progress_info_->start_scan(scan_, hmm_.iterations(), hmm_.connectivity(), blocks_.size(), epsilon);
        for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
            expectation_step(blocks_[b], epsilon);
            maximization_step();
            ++iterations_;
        }
        if (progress_info_) progress_info_->complete_scan(log_likelihood_, log_likelihood_change());
        if (epsilon < params_.epsilon_batch && static_cast<int>(blocks_.size()) > 1) {
            setup_blocks(true);
            epsilon = 1.0f;
        }
    }

    LOG(DEBUG) << hmm_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::expectation_step(const data_vector& block, float epsilon)
{
    // Run forward and backward algorithm on each subject in current block
    for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi) {
        ForwardBackwardMatrices fbm((*bi)->length(), hmm_.num_states());
        forward_backward_algorithm(hmm_, **bi, params_, fbm);

        if (progress_info_) progress_info_->progress(hmm_.num_states() * (**bi).length());
        add_contribution_to_priors(fbm);
        add_contribution_to_transitions(fbm);
        add_contribution_to_emissions(fbm, **bi);

        log_likelihood_ += fbm.log_likelihood / num_eff_cols_;
    }

    update_sufficient_statistics(epsilon);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::maximization_step()
{
    const int num_states    = hmm_.num_states();
    const int num_cols      = hmm_.num_cols();
    const int alphabet_size = hmm_[0].alphabet_size();

    // Calculate normalization factor for priors
    float sum = 0.0f;
    for (int k = 0; k < num_states; ++k) sum += profile_stats_[k]->prior();
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_[k];

        hmm_[k].set_prior(p_k.prior() * fac);
        ContextProfile<Alphabet_T> tmp(p_k);
        if (normalize(tmp)) {  // don't update profiles that did'n get any evidence
            tmp.transform_to_logspace();
            for (int i = 0; i < num_cols; ++i)
                for (int a = 0; a < alphabet_size; ++a)
                    hmm_[k][i][a] = tmp[i][a];
        }
    }

    // Calculate and assign new transition probabilities
    hmm_.clear_transitions();
    for (int k = 0; k < num_states; ++k) {
        sum = 0.0f;
        for (int l = 0; l < num_states; ++l) {
            if (transition_stats_.test(k,l)) {
                const float a_kl = transition_stats_[k][l] + params_.transition_pseudocounts - 1.0f;
                if (a_kl > 0.0f) {
                    transition_stats_[k][l] = a_kl;
                    sum += a_kl;
                } else {
                    transition_stats_.erase(k,l);
                }
            }
        }
        if (sum != 0.0f) {
            fac = 1.0f / sum;
            for (int l = 0; l < num_states; ++l) {
                if (transition_stats_.test(k,l)) {
                    hmm_(k,l) = transition_stats_[k][l] * fac;
                    LOG(DEBUG2) << strprintf("tr[%i][%i]=%-8.5f", k, l, static_cast<float>(hmm_(k,l)));
                }
            }
        }
    }

    // Increment iteration counter in HMM
    ++hmm_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_priors(const ForwardBackwardMatrices& m)
{
    const int num_states = hmm_.num_states();
    for (int k = 0; k < num_states; ++k)
        profile_stats_block_[k]->set_prior(profile_stats_block_[k]->prior() + m.f[0][k] * m.b[0][k]);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_transitions(const ForwardBackwardMatrices& m)
{
    const int slen = m.f.num_rows();

    for (const_transition_iterator ti = hmm_.transitions_begin(); ti != hmm_.transitions_end(); ++ti) {
        double w_kl = 0.0;
        for (int i = 0; i < slen-1; ++i) {
            w_kl += m.f[i][ti->from] * m.b[i+1][ti->to] * ti->probability * m.e[i+1][ti->to] / m.s[i+1];
        }

        if (w_kl != 0.0) {
            if (!transition_stats_block_.test(ti->from, ti->to)) transition_stats_block_[ti->from][ti->to] = 0.0f;
            transition_stats_block_[ti->from][ti->to] = transition_stats_block_[ti->from][ti->to] + w_kl;
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                                                                    const CountsProfile<Alphabet_T>& c)
{
    const int slen       = m.f.num_rows();
    const int num_states = hmm_.num_states();

    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];
        const int ci = p_k.center();

        for (int i = 0; i < slen; ++i) {
            const int beg = std::max(0, i - ci);
            const int end = std::min(c.num_cols() - 1, i + ci);

            for(int h = beg; h <= end; ++h) {
                const int j = h - i + ci;
                const int alphabet_size = p_k.alphabet_size();
                for (int a = 0; a < alphabet_size; ++a) {
                    p_k[j][a] += c[h][a] * m.f[i][k] * m.b[i][k];
                }
            }
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                                                                    const Sequence<Alphabet_T>& s)
{
    const int slen       = m.f.num_rows();
    const int num_states = hmm_.num_states();

    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];
        const int ci = p_k.center();

        for (int i = 0; i < slen; ++i) {
            const int beg = std::max(0, i - ci);
            const int end = std::min(s.length() - 1, i + ci);

            for(int h = beg; h <= end; ++h) {
                const int j = h - i + ci;
                p_k[j][s[h]] += m.f[i][k] * m.b[i][k];
            }
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::update_sufficient_statistics(float epsilon)
{
    const float gamma = 1.0f - epsilon;

    // Update transition statistics
    for (const_transition_iterator ti = hmm_.transitions_begin(); ti != hmm_.transitions_end(); ++ti) {
        if (transition_stats_block_.test(ti->from, ti->to)) {
            if (!transition_stats_.test(ti->from, ti->to)) transition_stats_[ti->from][ti->to] = 0.0f;
            transition_stats_[ti->from][ti->to] =
                gamma * transition_stats_[ti->from][ti->to] + transition_stats_block_[ti->from][ti->to];
        }
    }
    transition_stats_block_.clear();

    // Update priors and emissions statistics
    const int num_states = hmm_.num_states();
    const int num_cols   = hmm_.num_cols();
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_block = *profile_stats_block_[k];
        ContextProfile<Alphabet_T>& p       = *profile_stats_[k];
        const int alphabet_size = p.alphabet_size();

        p.set_prior(p.prior() * gamma + p_block.prior());
        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p[j][a] = gamma * p[j][a]  + p_block[j][a];
            }
        }
        reset(p_block, 0.0f);
        p_block.set_prior(0.0f);
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::init()
{
    // Create profiles for global and block-level sufficient statistics
    for (int k = 0; k < hmm_.num_states(); ++k) {
        profile_stats_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, hmm_.num_cols())));
        profile_stats_block_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, hmm_.num_cols())));
    }

    // Compute total number of data columns for progress bar
    int num_cols = 0;
    for (typename data_vector::const_iterator di = data_.begin(); di != data_.end(); ++di)
        num_cols += (**di).length();
    if (progress_info_) progress_info_->init(hmm_.num_states() * num_cols);

    // Set number of effective columsn for log-likelihood calculation
    ProfileMatcher<Alphabet_T> matcher;
    if (!params_.ignore_context && hmm_.num_cols() > 1)
        matcher.set_weights(hmm_.num_cols(), params_.weight_center, params_.weight_decay);
    num_eff_cols_ = matcher.num_eff_cols() * num_cols;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::setup_blocks(bool batch_mode)
{
    if (batch_mode) {
        if (static_cast<int>(blocks_.size()) == 1) return;  // already in batch mode
        blocks_.clear();
        blocks_.push_back(data_);

    } else {
        const int num_blocks = params_.num_blocks == 0 ? iround(pow(data_.size(), 3.0/8.0)) : params_.num_blocks;
        const int block_size = iround(data_.size() / num_blocks);

        for (int b = 0; b < num_blocks; ++b) {
            data_vector block;
            if (b == num_blocks - 1)  // last block may differ in block size
                for (int n = b * block_size; n < static_cast<int>(data_.size()); ++n)
                    block.push_back(data_[n]);
            else
                for (int n = b * block_size; n < (b+1) * block_size; ++n)
                    block.push_back(data_[n]);
            blocks_.push_back(block);
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline float BaumWelchTraining<Alphabet_T, Subject_T>::log_likelihood_change()
{
    return log_likelihood_ - log_likelihood_prev_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline bool BaumWelchTraining<Alphabet_T, Subject_T>::terminate()
{
    if (scan_ < params_.min_scans)
        return false;
    else if (scan_ >= params_.max_scans)
        return true;
    else if (params_.max_connectivity == 0)
        return fabs(log_likelihood_change()) <= params_.log_likelihood_change;
    else
        return fabs(log_likelihood_change()) <= params_.log_likelihood_change && hmm_.connectivity() <= params_.max_connectivity;
}

}  // cs

#endif
