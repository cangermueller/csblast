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
#include "hmm.h"
#include "log.h"
#include "context_profile.h"
#include "profile_matcher.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

template< class Alphabet_T>
class TrainingProgressInfo
{
  public:
    TrainingProgressInfo(const HMM<Alphabet_T>& hmm, std::ostream& out = std::cout)
            : hmm_(hmm),
              out_(out),
              total_(0),
              progress_(0),
              bar_(0)
    {
        out_ << strprintf("%-4s %-37s %10s %8s %12s\n", "Iter", "Progress", "log(L)", "+/-", "Connectivity");
        out_ << std::string(TABLE_WIDTH, '-') << std::endl;
    }

    void init(int total)
    {
        total_    = total;
        progress_ = 0;
        bar_      = 0;

        out_ << strprintf("%-3i  [", hmm_.iterations());
        out_.flush();
    }

    void increment(int incr)
    {
        progress_ += incr;
        const int bar_incr = round(static_cast<float>(progress_) / total_ * PRG_BAR_WIDTH) - bar_;
        std::string prg_str(bar_incr, '=');

        out_ << prg_str;
        if (progress_ == total_) out_ << "] 100% ";
        out_.flush();

        bar_ += bar_incr;
    }

    void print_stats(float log_likelihood, float delta = 1.0f, bool show_delta = true)
    {
        if (show_delta) {
            out_ << strprintf("%10.0f %+8.5f %12.2f\n", log_likelihood, delta, hmm_.connectivity());
            LOG(DEBUG) << strprintf("iter=%-3i log(L)=%-9.2e change=%-+8.5f connectivity=%-7.2f",
                                    hmm_.iterations(), log_likelihood, delta, hmm_.connectivity());
        } else {
            out_ << strprintf("%10.0f %8s %12.2f\n", log_likelihood, "", hmm_.connectivity());
            LOG(DEBUG) << strprintf("iter=%-3i log(L)=%-9.2e connectivity=%-7.2f",
                                    hmm_.iterations(), log_likelihood, hmm_.connectivity());
        }
        out_.flush();
    }

  private:
    static const int TABLE_WIDTH   = 75;
    static const int PRG_BAR_WIDTH = 30;

    // The HMM being trained.
    const HMM<Alphabet_T>& hmm_;
    // Output stream.
    std::ostream& out_;
    // Time complexity of iteration: O(NKL)
    int total_;
    // Progress so far.
    int progress_;
    // With of bar printed so far.
    int bar_;
};



struct BaumWelchParams : public ForwardBackwardParams
{
    BaumWelchParams()
            : ForwardBackwardParams(),
              min_scans(3),
              max_scans(100),
              max_connectivity(0),
              log_likelihood_threshold(0.001f),
              transition_pseudocounts(1.0f),
              num_blocks(1),
              epsilon_null(0.05f),
              beta(1.0f)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_scans(params.min_scans),
              max_scans(params.max_scans),
              max_connectivity(params.max_connectivity),
              log_likelihood_threshold(params.log_likelihood_threshold),
              transition_pseudocounts(params.transition_pseudocounts),
              num_blocks(params.num_blocks),
              epsilon_null(params.epsilon_null),
              beta(params.beta)
    { }

    // Minimal number of data scans.
    int min_scans;
    // Maximum number of data scans.
    int max_scans;
    // Maximum average connectivity for convergence.
    int max_connectivity;
    // Log-likelihood change per column for convergence.
    float log_likelihood_threshold;
    // Pseudocounts added to transitions (values below 1 enforce sparsity).
    float transition_pseudocounts;
    // Number of blocks into which the training data are divided (0 stands for automatics determination by B=N^(3/8)).
    int num_blocks;
    // Initial value for learning rate epsilon (1-epsilon is the preserved fraction of recently visited training profiles).
    float epsilon_null;
    // Paramter governing exponential decay of epsilon per scan.
    float beta;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchTraining
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
    typedef typename std::vector<data_vector> blocks_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;
    typedef typename HMM<Alphabet_T>::const_transition_iterator const_transition_iterator;

    // Initializes a new training object.
    BaumWelchTraining(const BaumWelchParams& params,
                      const data_vector& data,
                      HMM<Alphabet_T>& hmm,
                      TrainingProgressInfo<Alphabet_T>* progress_info = NULL);

    virtual ~BaumWelchTraining()
    { }

    // Trains the HMM with the data provided until one of the termination criterions is fullfilled.
    void run();

  private:
    // Prepares all members for HMM training.
    void init();
    // Fills the blocks vector with training data.
    void setup_blocks();
    // Runs forward backward algorithm on provided data.
    void run_forward_backward(const data_vector& block);
    // Adds the contribution of a subject's forward-backward matrices to prio probabilities of states.
    void add_contribution_to_priors(const ForwardBackwardMatrices& m);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Updates global sufficient statistics with sufficient statistics calculated on current block.
    void update_sufficient_statistics();
    // Calculates new HMM parameters by Maxmimum-Likelihood estimation.
    void assign_new_parameters();
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
    TrainingProgressInfo<Alphabet_T>* progress_info_;
    // Global expected sufficient statistics for transitions
    sparse_matrix<float> transition_stats_;
    // Global expeted sufficient statistics for emissions and state priors
    profiles_vector profile_stats_;
    // Expected sufficient statistics for transitions based on current block
    sparse_matrix<float> transition_stats_block_;
    // Expeted sufficient statistics for emissions and state priors based on current block
    profiles_vector profile_stats_block_;
    // Incremental likelihood in current iteration
    double log_likelihood_;
    // Complete likelihood after last complete scan of training data
    double log_likelihood_prev_;
    // Number of traning iterations performed so far.
    int iterations_;
    // Number of complete data scans performed so far.
    int scans_;
    // Total number of data subjects visited so far.
    int num_trained_;
    // Current learning rate epsilon (1-epsilon is the preserved fraction of recently visited training profiles).
    float epsilon_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline BaumWelchTraining<Alphabet_T, Subject_T>::BaumWelchTraining(const BaumWelchParams& params,
                                                                   const data_vector& data,
                                                                   HMM<Alphabet_T>& hmm,
                                                                   TrainingProgressInfo<Alphabet_T>* progress_info)
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
          scans_(0),
          num_trained_(0),
          epsilon_(params.epsilon_null)
{
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::init()
{
    for (int k = 0; k < hmm_.num_states(); ++k) {
        profile_stats_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, hmm_.num_cols())));
        profile_stats_block_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, hmm_.num_cols())));
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::setup_blocks()
{
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

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::run()
{
    LOG(DEBUG) << "Running Baum-Welch training on ...";
    LOG(DEBUG) << hmm_;

    // Calculate log-likelihood baseline by batch EM without blocks
    run_forward_backward(data_);
    update_sufficient_statistics();
    if (progress_info_) progress_info_->print_stats(log_likelihood_, log_likelihood_change(), false);

    setup_blocks();
    do {
        for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
            // Advance iteration counters and likelihoods
            ++hmm_;
            ++iterations_;
            log_likelihood_prev_ = log_likelihood_;
            log_likelihood_ = 0.0;

            assign_new_parameters();
            run_forward_backward(blocks_[b]);
            update_sufficient_statistics();

            if (progress_info_) progress_info_->print_stats(log_likelihood_, log_likelihood_change());
        }
        ++scans_;
    } while (!terminate());

    LOG(DEBUG) << hmm_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline bool BaumWelchTraining<Alphabet_T, Subject_T>::terminate()
{
    if (scans_ < params_.min_scans)
        return false;
    else if (scans_ >= params_.max_scans)
        return true;
    else if (params_.max_connectivity == 0)
        return fabs(log_likelihood_change()) <= params_.log_likelihood_threshold;
    else
        return fabs(log_likelihood_change()) <= params_.log_likelihood_threshold && hmm_.connectivity() <= params_.max_connectivity;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void BaumWelchTraining<Alphabet_T, Subject_T>::run_forward_backward(const data_vector& block)
{
    // Compute total number of columns in current block for progress bar
    if (progress_info_) {
        int num_cols_block = 0;
        for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi)
            num_cols_block += (**bi).length();
        progress_info_->init(hmm_.num_states() * num_cols_block);
    }

    // Run forward and backward algorithm on each subject in current block
    for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi) {
        ForwardBackwardMatrices fbm((*bi)->length(), hmm_.num_states());
        forward_backward_algorithm(hmm_, **bi, params_, fbm);
        if (progress_info_) progress_info_->increment(hmm_.num_states() * (**bi).length());

        add_contribution_to_priors(fbm);
        add_contribution_to_transitions(fbm);
        add_contribution_to_emissions(fbm, **bi);

        log_likelihood_ += fbm.log_likelihood;
    }
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
inline void BaumWelchTraining<Alphabet_T, Subject_T>::update_sufficient_statistics()
{
    // Update transition statistics
    for (const_transition_iterator ti = hmm_.transitions_begin(); ti != hmm_.transitions_end(); ++ti) {
        if (transition_stats_block_.test(ti->from, ti->to)) {
            if (!transition_stats_.test(ti->from, ti->to)) transition_stats_[ti->from][ti->to] = 0.0f;
            transition_stats_[ti->from][ti->to] =
                (1.0f - epsilon_) * transition_stats_[ti->from][ti->to] + transition_stats_block_[ti->from][ti->to];
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

        p.set_prior(p.prior() * (1.0f - epsilon_) + p_block.prior());
        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p[j][a] = (1.0f - epsilon_) * p[j][a]  + p_block[j][a];
            }
        }
        reset(p_block, 0.0f);
        p_block.set_prior(0.0f);
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::assign_new_parameters()
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
                if (transition_stats_[k][l] > 1.0f - params_.transition_pseudocounts) {
                    transition_stats_[k][l] = transition_stats_[k][l] + params_.transition_pseudocounts - 1.0f;
                    sum += transition_stats_[k][l];
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
    LOG(INFO) << hmm_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline float BaumWelchTraining<Alphabet_T, Subject_T>::log_likelihood_change()
{
    int data_cols = 0;
    for (typename data_vector::const_iterator di = data_.begin(); di != data_.end(); ++di)
        data_cols += (**di).length();
    return (log_likelihood_ - log_likelihood_prev_) / data_cols;
}

}  // cs

#endif
