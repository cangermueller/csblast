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
#include <valarray>
#include <vector>

#include "forward_backward_algorithm.h"  // FIXME: remove dependence on ForwardBackward params!
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
        out_ << std::string(76, '-') << std::endl;
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
    ClusteringParams()
            : ForwardBackwardParams(),
              min_scans(50),
              max_scans(500),
              log_likelihood_change(2e-4f),
              num_blocks(0),
              epsilon_null(0.5f),
              beta(0.2f),
              epsilon_batch(0.05f)
    { }

    ClusteringParams(const ClusteringParams& params)
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

    // Initializes a new clustering object.
    Clustering(const ClusteringParams& params,
               const data_vector& data,
               ProfileLibrary<Alphabet_T>& lib,
               ClusteringProgressTable* progress_table = NULL);

    virtual ~Clustering()
    { }

    // Optimize the profile library by EM clustering.
    void run();

  private:
    // Evaluates the responsibilities using the current parameter values.
    void expectation_step(const data_vector& block, float epsilon);
    // Reestimate teh parameters using the current responsibilities.
    void maximization_step();
    // Adds the contribution of the responsibilities for a subject to sufficient statistics for priors.
    void add_contribution_to_priors(const std::valarray<float>& p_zn);
    // Adds the contribution of the responsibilities for a counts profile to sufficient statistics for emissions.
    void add_contribution_to_emissions(const std::valarray<float>& p_zn, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of the responsibilities for a sequence to sufficient statistics for emissions.
    void add_contribution_to_emissions(const std::valarray<float>& p_zn, const Sequence<Alphabet_T>& s);
    // Updates global sufficient statistics with sufficient statistics calculated on current block.
    void update_sufficient_statistics(float epsilon);
    // Fills the blocks vector with training data.
    void setup_blocks(bool batch_mode = false);
    // Prepares all members for clustering.
    void init();
    // Calculates the change between current log-likelihood and the log-likelihood in previous scan.
    float log_likelihood_change();
    // Returns true if any termination condition is fullfilled.
    bool terminate();

    // Parameter wrapper for clustering.
    const ClusteringParams& params_;
    // Training data (either sequences or counts profiles)
    const data_vector& data_;
    // Blocks of training data.
    blocks_vector blocks_;
    // Profile library to be optimized.
    ProfileLibrary<Alphabet_T>& lib_;
    // Global expected sufficient statistics for emissions and state priors.
    profiles_vector profile_stats_;
    // Expeted sufficient statistics for emissions and state priors based on current block.
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
    // Profile matcher for calculation of emission probabilities.
    ProfileMatcher<Alphabet_T> profile_matcher_;
    // Progress info object for output of progress table.
    ClusteringProgressTable* progress_table_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline Clustering<Alphabet_T, Subject_T>::Clustering(const ClusteringParams& params,
                                                     const data_vector& data,
                                                     ProfileLibrary<Alphabet_T>& lib,
                                                     ClusteringProgressTable* progress_table)
        : params_(params),
          data_(data),
          blocks_(),
          lib_(lib),
          profile_stats_(),
          profile_stats_block_(),
          log_likelihood_(0.0f),
          log_likelihood_prev_(0.0f),
          iterations_(0),
          scan_(1),
          num_eff_cols_(0),
          profile_matcher_(),
          progress_table_(progress_table)
{
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::run()
{
    LOG(DEBUG) << "Running expecation-maximization clustering on ...";
    LOG(DEBUG) << lib_;

    if (progress_table_) progress_table_->print_header();

    // Calculate log-likelihood baseline by batch EM without blocks
    float epsilon = 1.0f;
    if (progress_table_) progress_table_->print_row_begin(scan_, lib_.iterations(), 1, epsilon);
    expectation_step(data_, epsilon);
    maximization_step();
    ++iterations_;
    if (progress_table_) progress_table_->print_row_end(log_likelihood_);

    // Perform E-step and M-step on each training block until convergence
    setup_blocks();
    while (!terminate()) {
        if (static_cast<int>(blocks_.size()) > 1) epsilon = params_.epsilon_null * exp(-params_.beta * (scan_ - 1));
        log_likelihood_prev_ = log_likelihood_;
        log_likelihood_      = 0.0;
        ++scan_;

        if (progress_table_) progress_table_->print_row_begin(scan_, lib_.iterations(), blocks_.size(), epsilon);
        for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
            expectation_step(blocks_[b], epsilon);
            maximization_step();
            ++iterations_;
        }
        if (progress_table_) progress_table_->print_row_end(log_likelihood_, log_likelihood_change());
        if (epsilon < params_.epsilon_batch && static_cast<int>(blocks_.size()) > 1) {
            setup_blocks(true);
            epsilon = 1.0f;
        }
    }

    LOG(DEBUG) << lib_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::expectation_step(const data_vector& block, float epsilon)
{
    const int num_profiles = lib_.num_profiles();
    std::valarray<float> p_zn(0.0f, lib_.num_profiles());

    // Given each training window compute posterior probabilities p_zn[k] of profile k
    for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi) {
        float sum = 0.0f;
        for (int k = 0; k < num_profiles; ++k) {
            p_zn[k] = lib_[k].prior() * profile_matcher_(lib_[k], **bi, lib_[k].center());
            sum += p_zn[k];
        }
        p_zn /= sum;
        add_contribution_to_priors(p_zn);
        add_contribution_to_emissions(p_zn, **bi);
        log_likelihood_ += log(sum) / num_eff_cols_;

        if (progress_table_) progress_table_->print_progress(lib_.num_profiles());
    }

    update_sufficient_statistics(epsilon);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::maximization_step()
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    float sum = 0.0f;
    for (int k = 0; k < num_profiles; ++k) sum += profile_stats_[k]->prior();
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_[k];

        lib_[k].set_prior(p_k.prior() * fac);
        ContextProfile<Alphabet_T> tmp(p_k);
        if (normalize(tmp)) {  // don't update profiles that did'n get any evidence
            tmp.transform_to_logspace();
            for (int i = 0; i < num_cols; ++i)
                for (int a = 0; a < alphabet_size; ++a)
                    lib_[k][i][a] = tmp[i][a];
        }
    }

    // Increment iteration counter in HMM
    ++lib_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_priors(const std::valarray<float>& p_zn)
{
    const int num_profiles = lib_.num_profiles();
    for (int k = 0; k < num_profiles; ++k)
        profile_stats_block_[k]->set_prior(profile_stats_block_[k]->prior() + p_zn[k]);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_emissions(const std::valarray<float>& p_zn,
                                                                             const CountsProfile<Alphabet_T>& c)
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];

        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p_k[j][a] += c[j][a] * p_zn[k];
            }
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_emissions(const std::valarray<float>& p_zn,
                                                                             const Sequence<Alphabet_T>& s)
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];

        for (int j = 0; j < num_cols; ++j) {
            p_k[j][s[j]] += p_zn[k];
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::update_sufficient_statistics(float epsilon)
{
    const float gamma       = 1.0f - epsilon;
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    // Update priors and emissions statistics
    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_block = *profile_stats_block_[k];
        ContextProfile<Alphabet_T>& p       = *profile_stats_[k];

        p.set_prior(p.prior() * gamma + p_block.prior());
        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p[j][a] = gamma * p[j][a]  + p_block[j][a];
            }
        }
        reset(p_block);
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::init()
{
    // Create profiles for global and block-level sufficient statistics
    for (int k = 0; k < lib_.num_profiles(); ++k) {
        profile_stats_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, lib_.num_cols())));
        profile_stats_block_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, lib_.num_cols())));
    }

    // Initialize total amount of work per scan
    if (progress_table_) progress_table_->init(lib_.num_profiles() * data_.size());

    // Set number of effective columns for log-likelihood calculation
    if (!params_.ignore_context && lib_.num_cols() > 1)
        profile_matcher_.set_weights(lib_.num_cols(), params_.weight_center, params_.weight_decay);
    num_eff_cols_ = profile_matcher_.num_eff_cols() * data_.size();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::setup_blocks(bool batch_mode)
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
inline float Clustering<Alphabet_T, Subject_T>::log_likelihood_change()
{
    return log_likelihood_ - log_likelihood_prev_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline bool Clustering<Alphabet_T, Subject_T>::terminate()
{
    if (scan_ < params_.min_scans)
        return false;
    else if (scan_ >= params_.max_scans)
        return true;
    else
        return fabs(log_likelihood_change()) <= params_.log_likelihood_change;
}

}  // cs

#endif
