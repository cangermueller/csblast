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
              delta_log_likelihood_threshold_(1e-4)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_iterations_(params.min_iterations_),
              max_iterations_(params.max_iterations_),
              delta_log_likelihood_threshold_(params.delta_log_likelihood_threshold_)
    { }

    float min_iterations() const { return min_iterations_; }
    float max_iterations() const { return max_iterations_; }
    float delta_log_likelihood_threshold() const { return delta_log_likelihood_threshold_; }
    void min_iterations(float min_iter) { min_iterations_ = min_iter; }
    void max_iterations(float max_iter) { max_iterations_ = max_iter; }
    void delta_log_likelihoodthreshold_(float threshold) { delta_log_likelihood_threshold_ = threshold; }

  private:
    float min_iterations_;
    float max_iterations_;
    float delta_log_likelihood_threshold_;
};



class TrainingProgressInfo
{
  public:
    TrainingProgressInfo(std::ostream& out = std::cout)
            : out_(out),
              total_(0),
              progress_(0),
              bar_(0),
              iter_(0)
    {
        out_ << strprintf("%-4s %-47s  %7s %8s\n", "Iter", "Progress", "-log(L)", "+/-\%");
        out_ << std::string(TABLE_WIDTH, '-') << std::endl;
    }

    void init(int iter, int total)
    {
        iter_     = iter;
        total_    = total;
        progress_ = 0;
        bar_      = 0;

        out_ << strprintf("%-3i  [", iter);
        out_.flush();
    }

    void increment(int incr)
    {
        progress_ += incr;
        const int bar_incr = round(static_cast<float>(progress_) / total_ * PRG_BAR_WIDTH) - bar_;
        std::string prg_str(bar_incr, '=');

        out_ << prg_str;
        if (progress_ == total_) out_ << "] 100%  ";
        out_.flush();

        bar_ += bar_incr;
    }

    void print_stats(float log_likelihood, float delta)
    {
        if (iter_ == 0) {
            out_ << strprintf("%7.0f\n", log_likelihood);
            LOG(DEBUG) << strprintf("%-3i  %7.0f", iter_, log_likelihood);
        } else {
            out_ << strprintf("%7.0f %+7.3f", log_likelihood, delta * 100.0f) << "%\n";
            LOG(DEBUG) << strprintf("%-3i  %7.0f %+7.3f", iter_, log_likelihood, delta * 100.0f) << "%";
        }
        out_.flush();
    }

  private:
    static const int TABLE_WIDTH   = 70;
    static const int PRG_BAR_WIDTH = 40;

    // Output stream.
    std::ostream& out_;
    // Time complexity of iteration: O(NKL)
    int total_;
    // Progress so far.
    int progress_;
    // With of bar printed so far.
    int bar_;
    // .Current iteration
    int iter_;
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
    float run(HMM<Alphabet_T>& hmm, const data_vector& data, TrainingProgressInfo* prg_info = NULL);

  private:
    // Prepares the stage for a new training run.
    void setup(int num_states, int num_cols);
    // Runs forward backward algorithm on provided data.
    void run_forward_backward(HMM<Alphabet_T>& hmm, const data_vector& data, TrainingProgressInfo* prg_info = NULL);
    // Adds the contribution of a subject's forward-backward matrices to prio probabilities of states.
    void add_contribution_to_priors(const ForwardBackwardMatrices& m);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m, const HMM<Alphabet_T>& hmm);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Calculates new HMM parameters by Maxmimum-Likelihood estimation.
    void calculate_and_apply_new_parameters(HMM<Alphabet_T>& hmm);

    // Calculates the the percent difference between the current log-likelihood and the previous-iteration log-likelihood.
    float delta_log_likelihood()
    {
        return log_likelihood_prev_ == 0.0f ? 1.0f : (log_likelihood_ - log_likelihood_prev_) / log_likelihood_prev_;
    }

    // Instance of forward-backward algorithm to compute exptected number of transitions and emissions.
    ForwardBackwardAlgorithm<Alphabet_T, Subject_T>* fb_;
    // Transition counts
    sparse_matrix<float> transitions_;
    // Profiles with emission counts
    profiles_vector profiles_;
    // Likelihood of current iteration
    double log_likelihood_;
    // Likelihood of previous iteration
    double log_likelihood_prev_;
    // Number of traning iterations applied so far.
    int iter_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
float BaumWelchTraining<Alphabet_T, Subject_T>::run(HMM<Alphabet_T>& hmm,
                                                    const data_vector& data,
                                                    TrainingProgressInfo* prg_info)
{
    LOG(DEBUG) << "Running Baum-Welch training on ...";
    LOG(DEBUG) << hmm;
    setup(hmm.num_states(), hmm[1].num_cols());

    // Calculate log-likelihood baseline
    run_forward_backward(hmm, data, prg_info);
    if (prg_info) prg_info->print_stats(log_likelihood_, delta_log_likelihood());

    do {
        calculate_and_apply_new_parameters(hmm);
        run_forward_backward(hmm, data, prg_info);

        if (prg_info) prg_info->print_stats(log_likelihood_, delta_log_likelihood());

    } while(iter_ < min_iterations() || iter_ < max_iterations() &&
            fabs(delta_log_likelihood()) > delta_log_likelihood_threshold());

    LOG(INFO) << hmm;
    return log_likelihood_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::run_forward_backward(HMM<Alphabet_T>& hmm,
                                                                    const data_vector& data,
                                                                    TrainingProgressInfo* prg_info)
{
    if (prg_info) {
        int total = 0;
        for (typename data_vector::const_iterator di = data.begin(); di != data.end(); ++di)
            total += hmm.num_states() * (**di).length();
        prg_info->init(iter_, total);
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
    iter_                = 0;
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
    const int slen       = m.f.num_rows();
    for (const_transition_iterator ti = hmm.transitions_begin(); ti != hmm.transitions_end(); ++ti) {
        if (!transitions_.test(ti->from, ti->to)) transitions_[ti->from][ti->to] = 0.0f;

        double prob_a_kl = 0.0f;
        for (int i = 0; i < slen-1; ++i) {
            prob_a_kl += m.s[i+1] * m.f[i][ti->from] * m.b[i+1][ti->to] * ti->probability * m.e[i+1][ti->to];
        }
        transitions_[ti->from][ti->to] = transitions_[ti->from][ti->to] + prob_a_kl;
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
void BaumWelchTraining<Alphabet_T, Subject_T>::calculate_and_apply_new_parameters(HMM<Alphabet_T>& hmm)
{
    // Prepare for next iteration
    ++iter_;
    log_likelihood_prev_ = log_likelihood_;
    log_likelihood_ = 0.0;

    // Calculate normalization factor for priors and normalize profiles
    const int num_states = hmm.num_states();
    float sum = 0.0f;
    for (int k = 0; k < num_states; ++k) {
        sum += profiles_[k]->prior();
        normalize(*profiles_[k]);
    }
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_states; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];

        hmm[k].set_prior(p_k.prior() * fac);

        p_k.transform_to_logspace();
        const int num_cols      = p_k.num_cols();
        const int alphabet_size = p_k.alphabet_size();
        for (int i = 0; i < num_cols; ++i)
            for (int a = 0; a < alphabet_size; ++a)
                hmm[k][i][a] = p_k[i][a];
        p_k.transform_to_linspace();

        reset(p_k, 0.0f);
        p_k.set_prior(0.0f);
    }

    // Calculate and assign new transition probabilities
    for (int k = 0; k < num_states; ++k) {
        sum = 0.0f;
        for (int l = 0; l < num_states; ++l)
            if (transitions_.test(k,l)) sum += transitions_[k][l];

        if (sum != 0.0f) {
            fac = 1.0f / sum;
            for (int l = 0; l < num_states; ++l)
                if (transitions_[k][l] == 0.0f) {
                    transitions_.erase(k,l);
                } else if (transitions_.test(k,l)) {
                    hmm(k,l) = transitions_[k][l] * fac;
                }
        } else {
            throw Exception("Unable calculate new HMM transition probabilities: state %i has no out-transitions!", k);
        }
    }
    transitions_.clear();
}

}  // cs

#endif
