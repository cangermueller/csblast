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

#include "context_profile-inl.h"
#include "emitter.h"
#include "expectation_maximization.h"
#include "forward_backward_algorithm.h"
#include "hmm.h"
#include "log.h"
#include "progress_table.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

// Forwrad declarations
template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchProgressTable;


struct BaumWelchParams : public EmissionParams, public ExpectationMaximizationParams
{
    BaumWelchParams()
            : EmissionParams(),
              ExpectationMaximizationParams(),
              transition_pseudocounts(1.0f),
              max_connectivity(0)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : EmissionParams(params),
              ExpectationMaximizationParams(params),
              transition_pseudocounts(params.transition_pseudocounts),
              max_connectivity(params.max_connectivity)
    { }

    // Pseudocounts added to transitions (values below 1 enforce sparsity).
    float transition_pseudocounts;
    // Maximum average connectivity for convergence.
    int max_connectivity;
};


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchTraining : public ExpectationMaximization<Alphabet_T, Subject_T>
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;
    typedef typename HMM<Alphabet_T>::const_transition_iterator const_transition_iterator;

    // Needed to access names in templatized base class
    using ExpectationMaximization<Alphabet_T, Subject_T>::log_likelihood_change;

    // Initializes a new training object without output.
    BaumWelchTraining(const BaumWelchParams& params,
                      const data_vector& data,
                      HMM<Alphabet_T>& hmm);
    // Initializes a new training object with output.
    BaumWelchTraining(const BaumWelchParams& params,
                      const data_vector& data,
                      HMM<Alphabet_T>& hmm,
                      std::ostream& out);

    virtual ~BaumWelchTraining();

  protected:
    // Needed to access names in templatized base class
    using ExpectationMaximization<Alphabet_T, Subject_T>::progress_table_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::log_likelihood_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::num_eff_cols_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::epsilon_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::data_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::scan_;

    // Runs forward backward algorithm on provided data.
    virtual void expectation_step(const data_vector& block);
    // Calculates and assigns new HMM parameters by Maxmimum-Likelihood estimation.
    virtual void maximization_step();
    // Prepares all members for HMM training.
    virtual void init();
    // Returns true if any termination condition is fullfilled.
    virtual bool terminate() const;
    // Returns parameter wrapper
    virtual const BaumWelchParams& params() const { return params_; }
    // Adds the contribution of a subject's forward-backward matrices to prio probabilities of states.
    void add_contribution_to_priors(const ForwardBackwardMatrices& m);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Updates global sufficient statistics with sufficient statistics calculated on current block.
    void update_sufficient_statistics();

    // Parameter wrapper for clustering.
    const BaumWelchParams& params_;
    // HMM to be trained
    HMM<Alphabet_T>& hmm_;
    // Profile matcher for calculation of emission probabilities.
    Emitter<Alphabet_T> emitter_;
    // Global expected sufficient statistics for transitions
    sparse_matrix<float> transition_stats_;
    // Global expeted sufficient statistics for emissions and state priors
    profiles_vector profile_stats_;
    // Expected sufficient statistics for transitions based on current block
    sparse_matrix<float> transition_stats_block_;
    // Expeted sufficient statistics for emissions and state priors based on current block
    profiles_vector profile_stats_block_;

    friend class BaumWelchProgressTable<Alphabet_T, Subject_T>;
};


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchProgressTable : public ProgressTable
{
  public:
    BaumWelchProgressTable(const BaumWelchTraining<Alphabet_T, Subject_T>* training,
                           std::ostream& out = std::cout,
                           int width = 30);

    virtual ~BaumWelchProgressTable() { }

    virtual void print_header();
    virtual void print_row_begin();
    virtual void print_row_end();

  protected:
    const BaumWelchTraining<Alphabet_T, Subject_T>* training_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
BaumWelchTraining<Alphabet_T, Subject_T>::BaumWelchTraining(const BaumWelchParams& params,
                                                            const data_vector& data,
                                                            HMM<Alphabet_T>& hmm)
        : ExpectationMaximization<Alphabet_T, Subject_T>(params, data),
          params_(params),
          hmm_(hmm),
          emitter_(hmm.num_cols(), params),
          transition_stats_(hmm.num_states(), hmm.num_states()),
          profile_stats_(),
          transition_stats_block_(hmm.num_states(), hmm.num_states()),
          profile_stats_block_()
{
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
BaumWelchTraining<Alphabet_T, Subject_T>::BaumWelchTraining(const BaumWelchParams& params,
                                                            const data_vector& data,
                                                            HMM<Alphabet_T>& hmm,
                                                            std::ostream& out)
        : ExpectationMaximization<Alphabet_T, Subject_T>(data),
          params_(params),
          hmm_(hmm),
          emitter_(hmm.num_cols(), params),
          transition_stats_(hmm.num_states(), hmm.num_states()),
          profile_stats_(),
          transition_stats_block_(hmm.num_states(), hmm.num_states()),
          profile_stats_block_()
{
    progress_table_ = new BaumWelchProgressTable<Alphabet_T, Subject_T>(this, out);
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
BaumWelchTraining<Alphabet_T, Subject_T>::~BaumWelchTraining()
{
    if (progress_table_) delete progress_table_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::expectation_step(const data_vector& block)
{
    // Run forward and backward algorithm on each subject in current block
    for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi) {
        ForwardBackwardMatrices fbm((*bi)->length(), hmm_.num_states());
        forward_backward_algorithm(hmm_, **bi, emitter_, &fbm);
        add_contribution_to_priors(fbm);
        add_contribution_to_transitions(fbm);
        add_contribution_to_emissions(fbm, **bi);
        log_likelihood_ += fbm.log_likelihood / num_eff_cols_;

        if (progress_table_) progress_table_->print_progress(hmm_.num_states() * (**bi).length());
    }

    update_sufficient_statistics();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::maximization_step()
{
    const int num_states    = hmm_.num_states();
    const int num_cols      = hmm_.num_cols();
    const int alphabet_size = hmm_.alphabet_size();

    // Calculate normalization factor for priors
    float sum = 0.0f;
    for (int k = 0; k < num_states; ++k) sum += profile_stats_[k]->prior();
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_[k];

        hmm_[k].set_prior(p_k.prior() * fac);
        ContextProfile<Alphabet_T> tmp(p_k);
        if (normalize(&tmp)) {  // don't update profiles that did'n get any evidence
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
                                                                                    const CountProfile<Alphabet_T>& c)
{
    const int slen       = c.length();
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
    const int slen       = s.length();
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
    const float gamma       = 1.0f - epsilon_;
    const int num_states    = hmm_.num_states();
    const int num_cols      = hmm_.num_cols();
    const int alphabet_size = hmm_.alphabet_size();

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
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_block = *profile_stats_block_[k];
        ContextProfile<Alphabet_T>& p       = *profile_stats_[k];

        p.set_prior(p.prior() * gamma + p_block.prior());
        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p[j][a] = gamma * p[j][a]  + p_block[j][a];
            }
        }
        reset(&p_block);
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

    // Compute total number of data columns
    int num_cols = 0;
    for (typename data_vector::const_iterator di = data_.begin(); di != data_.end(); ++di)
        num_cols += (**di).length();
    if (progress_table_) progress_table_->set_total_work(hmm_.num_states() * num_cols);

    // Set number of effective columsn for log-likelihood calculation
    num_eff_cols_ = emitter_.sum_weights() * num_cols;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline bool BaumWelchTraining<Alphabet_T, Subject_T>::terminate() const
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



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
BaumWelchProgressTable<Alphabet_T, Subject_T>::BaumWelchProgressTable(const BaumWelchTraining<Alphabet_T, Subject_T>* training,
                                                                      std::ostream& out,
                                                                      int width)
        : ProgressTable(out, width),
          training_(training)
{ }

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchProgressTable<Alphabet_T, Subject_T>::print_header()
{
    out_ << strprintf("%-4s %4s %6s %4s %7s  %-30s  %9s  %8s\n",
                      "Scan", "Itrs", "Conn", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
    out_ << std::string(82, '-') << std::endl;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchProgressTable<Alphabet_T, Subject_T>::print_row_begin()
{
    reset();
    out_ << strprintf("%-4i %4i %6.1f %4i %7.4f  ", training_->scan(), training_->iterations(),
                      training_->hmm_.connectivity(), training_->num_blocks(), training_->epsilon());
    out_.flush();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchProgressTable<Alphabet_T, Subject_T>::print_row_end()
{
    if (training_->scan() == 1)
        out_ << strprintf("  %9.5f\n", training_->log_likelihood());
    else
        out_ << strprintf("  %9.5f  %+8.5f\n", training_->log_likelihood(), training_->log_likelihood_change());
    out_.flush();
}

}  // cs

#endif
