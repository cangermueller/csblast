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
#include "util.h"

namespace cs
{

class BaumWelchParams : public ForwardBackwardParams
{
  public:
    BaumWelchParams()
            : ForwardBackwardParams(),
              min_iterations_(3),
              max_iterations_(100),
              max_connectivity_(0),
              log_likelihood_threshold_(0.001f),
              alpha_(1.0f)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_iterations_(params.min_iterations_),
              max_iterations_(params.max_iterations_),
              max_connectivity_(params.max_connectivity_),
              log_likelihood_threshold_(params.log_likelihood_threshold_),
              alpha_(params.alpha_)
    { }

    virtual ~BaumWelchParams()
    { }

    int min_iterations() const { return min_iterations_; }
    int max_iterations() const { return max_iterations_; }
    int max_connectivity() const { return max_connectivity_; }
    float log_likelihood_threshold() const { return log_likelihood_threshold_; }
    float alpha() const { return alpha_; }

    BaumWelchParams& min_iterations(int min_iter) { min_iterations_ = min_iter; return *this; }
    BaumWelchParams& max_iterations(int max_iter) { max_iterations_ = max_iter; return *this; }
    BaumWelchParams& max_connectivity(int max_connect) { max_connectivity_ = max_connect; return *this; }
    BaumWelchParams& log_likelihood_threshold(float t) { log_likelihood_threshold_ = t; return *this; }
    BaumWelchParams& alpha(float pc) { alpha_ = pc; return *this; }

  protected:
    int min_iterations_;
    int max_iterations_;
    int max_connectivity_;
    float log_likelihood_threshold_;
    float alpha_;
};


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



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchTraining : public BaumWelchParams
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;
    typedef typename HMM<Alphabet_T>::const_transition_iterator const_transition_iterator;

    BaumWelchTraining()
            : fb_(NULL),
              log_likelihood_(0.0f),
              log_likelihood_prev_(0.0f)
    { }

    BaumWelchTraining(const BaumWelchParams& params)
            : BaumWelchParams(params),
              fb_(NULL),
              log_likelihood_(0.0f),
              log_likelihood_prev_(0.0f)
    { }

    virtual ~BaumWelchTraining()
    {
        if (fb_) delete fb_;
    }

    // Trains the HMM with the data provided until one of the termination criterions is fullfilled.
    void run(HMM<Alphabet_T>& hmm,
             const data_vector& data,
             TrainingProgressInfo<Alphabet_T>* prg_info = NULL);

  private:
    // Prepares the stage for a new training run.
    void setup(int num_states, int num_cols);
    // Runs forward backward algorithm on provided data.
    void run_forward_backward(HMM<Alphabet_T>& hmm,
                              const data_vector& data,
                              TrainingProgressInfo<Alphabet_T>* prg_info = NULL);
    // Adds the contribution of a subject's forward-backward matrices to prio probabilities of states.
    void add_contribution_to_priors(const ForwardBackwardMatrices& m);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m, const HMM<Alphabet_T>& hmm);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Calculates new HMM parameters by Maxmimum-Likelihood estimation.
    void assign_new_parameters(HMM<Alphabet_T>& hmm);
    // Calculates the change between the current log-likelihood and the previous-iteration log-likelihood.
    float log_likelihood_change(const data_vector& data);

    // Instance of forward-backward algorithm to compute exptected number of transitions and emissions.
    ForwardBackwardAlgorithm<Alphabet_T, Subject_T>* fb_;
    // Transition counts
    sparse_matrix<float> transitions_;
    // Profiles with emission counts and prior counts
    profiles_vector profiles_;
    // Likelihood of current iteration
    double log_likelihood_;
    // Likelihood of previous iteration
    double log_likelihood_prev_;
    // Number of traning iterations applied so far.
    int iterations_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::run(HMM<Alphabet_T>& hmm,
                                                   const data_vector& data,
                                                   TrainingProgressInfo<Alphabet_T>* prg_info)
{
    LOG(DEBUG) << "Running Baum-Welch training on ...";
    LOG(DEBUG) << hmm;
    setup(hmm.num_states(), hmm.num_cols());

    // Calculate log-likelihood baseline
    run_forward_backward(hmm, data, prg_info);
    if (prg_info) prg_info->print_stats(log_likelihood_, log_likelihood_change(data), false);

    do {
        assign_new_parameters(hmm);
        run_forward_backward(hmm, data, prg_info);

        if (prg_info) prg_info->print_stats(log_likelihood_, log_likelihood_change(data));

    } while( iterations_ < min_iterations() || iterations_ < max_iterations() &&
             (fabs(log_likelihood_change(data)) > log_likelihood_threshold() ||
             max_connectivity() != 0 && hmm.connectivity() > max_connectivity()) );

    LOG(DEBUG) << hmm;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::run_forward_backward(HMM<Alphabet_T>& hmm,
                                                                    const data_vector& data,
                                                                    TrainingProgressInfo<Alphabet_T>* prg_info)
{
    if (prg_info) {
        int total = 0;
        for (typename data_vector::const_iterator di = data.begin(); di != data.end(); ++di)
            total += hmm.num_states() * (**di).length();
        prg_info->init(total);
    }

    for (typename data_vector::const_iterator di = data.begin(); di != data.end(); ++di) {
            shared_ptr<ForwardBackwardMatrices> fbm = fb_->run(hmm, **di);
            if (prg_info) prg_info->increment(hmm.num_states() * (**di).length());

            add_contribution_to_priors(*fbm);
            add_contribution_to_transitions(*fbm, hmm);
            add_contribution_to_emissions(*fbm, **di);

            log_likelihood_ += fbm->log_likelihood;
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::setup(int num_states, int num_cols)
{
    if (fb_) delete fb_;
    fb_ = new ForwardBackwardAlgorithm<Alphabet_T, Subject_T>(*this);

    transitions_.clear();
    transitions_.resize(num_states, num_states);

    profiles_.clear();
    for (int k = 0; k < num_states; ++k) {
        shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(k, num_cols));
        profiles_.push_back(profile_ptr);
    }

    log_likelihood_      = 0.0f;
    log_likelihood_prev_ = 0.0f;
    iterations_          = 0;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_priors(const ForwardBackwardMatrices& m)
{
    const int num_states = profiles_.size();
    for (int k = 0; k < num_states; ++k)
        profiles_[k]->set_prior(profiles_[k]->prior() + m.f[0][k] * m.b[0][k]);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_transitions(const ForwardBackwardMatrices& m,
                                                                               const HMM<Alphabet_T>& hmm)
{
    const int slen = m.f.num_rows();
    for (const_transition_iterator ti = hmm.transitions_begin(); ti != hmm.transitions_end(); ++ti) {
        double w_kl = 0.0;
        for (int i = 0; i < slen-1; ++i) {
            w_kl += m.f[i][ti->from] * m.b[i+1][ti->to] * ti->probability * m.e[i+1][ti->to] / m.s[i+1];
        }

        if (w_kl != 0.0) {
            if (!transitions_.test(ti->from, ti->to)) transitions_[ti->from][ti->to] = 0.0f;
            transitions_[ti->from][ti->to] = transitions_[ti->from][ti->to] + w_kl;
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                                                             const CountsProfile<Alphabet_T>& c)
{
    const int slen       = m.f.num_rows();
    const int num_states = transitions_.num_rows();
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];
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
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                                                             const Sequence<Alphabet_T>& s)
{
    const int slen       = m.f.num_rows();
    const int num_states = transitions_.num_rows();

    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];
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
void BaumWelchTraining<Alphabet_T, Subject_T>::assign_new_parameters(HMM<Alphabet_T>& hmm)
{
    const int num_states    = hmm.num_states();
    const int num_cols      = hmm.num_cols();
    const int alphabet_size = hmm[0].alphabet_size();

    // Advance iteration counters and likelihoods
    ++hmm;
    ++iterations_;
    log_likelihood_prev_ = log_likelihood_;
    log_likelihood_ = 0.0;

    // Calculate normalization factor for priors
    float sum = 0.0f;
    for (int k = 0; k < num_states; ++k) sum += profiles_[k]->prior();
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];

        hmm[k].set_prior(p_k.prior() * fac);
        if (normalize(p_k)) {  // don't update profiles that did'n get any evidence
            p_k.transform_to_logspace();
            for (int i = 0; i < num_cols; ++i)
                for (int a = 0; a < alphabet_size; ++a)
                    hmm[k][i][a] = p_k[i][a];
            p_k.transform_to_linspace();
        }
        reset(p_k, 0.0f);
        p_k.set_prior(0.0f);
    }

    // Calculate and assign new transition probabilities
    hmm.clear_transitions();
    for (int k = 0; k < num_states; ++k) {
        sum = 0.0f;
        for (int l = 0; l < num_states; ++l) {
            if (transitions_.test(k,l)) {
                if (transitions_[k][l] > 1.0f - alpha()) {
                    transitions_[k][l] = transitions_[k][l] - 1.0f + alpha();
                    sum += transitions_[k][l];
                } else {
                    transitions_.erase(k,l);
                }
            }
        }
        if (sum != 0.0f) {
            fac = 1.0f / sum;
            for (int l = 0; l < num_states; ++l) {
                if (transitions_.test(k,l)) {
                    hmm(k,l) = transitions_[k][l] * fac;
                    LOG(DEBUG2) << strprintf("tr[%i][%i]=%-8.5f", k, l, static_cast<float>(hmm(k,l)));
                }
            }
        }
    }
    LOG(INFO) << hmm;
    transitions_.clear();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
float BaumWelchTraining<Alphabet_T, Subject_T>::log_likelihood_change(const data_vector& data)
{
    int data_cols = 0;
    for (typename data_vector::const_iterator di = data.begin(); di != data.end(); ++di)
        data_cols += (**di).length();
    return (log_likelihood_ - log_likelihood_prev_) / data_cols;
}

}  // cs

#endif
