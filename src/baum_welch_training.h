#ifndef CS_BAUM_WELCH_TRAINING_H
#define CS_BAUM_WELCH_TRAINING_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of Baum-Welch training for HMMs.

#include <cmath>

#include <vector>

#include "forward_backward_algorithm.h"
#include "hmm.h"
#include "log.h"
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
              max_iterations_(10),
              delta_likelihood_threshold_(1e-4)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_iterations_(params.min_iterations_),
              max_iterations_(params.max_iterations_),
              delta_likelihood_threshold_(params.delta_likelihood_threshold_)
    { }

    float min_iterations() const { return min_iterations_; }
    float max_iterations() const { return max_iterations_; }
    float delta_likelihood_threshold() const { return delta_likelihood_threshold_; }
    void min_iterations(float min_iter) { min_iterations_ = min_iter; }
    void max_iterations(float max_iter) { max_iterations_ = max_iter; }
    void delta_likelihoodthreshold_(float threshold) { delta_likelihood_threshold_ = threshold; }

  private:
    float min_iterations_;
    float max_iterations_;
    float delta_likelihood_threshold_;
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
              likelihood_(0.0f),
              likelihood_prev_(0.0f)
    { }

    BaumWelchTraining(const BaumWelchParams& params)
            : BaumWelchParams(params),
              fb_(NULL),
              likelihood_(0.0f),
              likelihood_prev_(0.0f)
    { }

    virtual ~BaumWelchTraining()
    {
        if (fb_) delete fb_;
    }

    // Trains the HMM with the data provided until one of the termination criterions is fullfilled.
    float run(HMM<Alphabet_T>& hmm, const data_vector& data);

  private:
    // Prepares the stage for a new training run.
    void setup(int num_states, int num_cols);
    // Adds the contribution of a subject's forward-backward matrices to transition counts.
    void add_contribution_to_transitions(const ForwardBackwardMatrices& m, const HMM<Alphabet_T>& hmm);
    // Adds the contribution of a count profile's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of a sequence's forward-backward matrices to emission counts.
    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s);
    // Calculates new HMM parameters by Maxmimum-Likelihood estimation.
    void calculate_and_apply_new_parameters(HMM<Alphabet_T>& hmm);

    // Calculates the the percent difference between the current likelihood and the previous-iteration likelihood.
    float delta_likelihood()
    {
        return likelihood_prev_ == 0.0f ? 1.0 : fabs(likelihood_ - likelihood_prev_) / likelihood_prev_;
    }

    // Instance of forward-backward algorithm to compute exptected number of transitions and emissions.
    ForwardBackwardAlgorithm<Alphabet_T, Subject_T>* fb_;
    // Transition counts
    sparse_matrix<float> transitions_;
    // Profiles with emission counts
    profiles_vector profiles_;
    // Likelihood of current iteration
    double likelihood_;
    // Likelihood of previous iteration
    double likelihood_prev_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
float BaumWelchTraining<Alphabet_T, Subject_T>::run(HMM<Alphabet_T>& hmm, const data_vector& data)
{
    LOG(INFO) << "Running Baum-Welch training on ...";
    LOG(INFO) << hmm;
    setup(hmm.num_states(), hmm[1].num_cols());

    for (int i = 1; i <= min_iterations() || i <= max_iterations() && delta_likelihood() > delta_likelihood_threshold(); ++i) {
        likelihood_prev_ = likelihood_;
        likelihood_ = 0.0f;

        for (typename data_vector::const_iterator di = data.begin(); di != data.end(); ++di) {
            shared_ptr<ForwardBackwardMatrices> fbm = fb_->run(hmm, **di);
            add_contribution_to_transitions(*fbm, hmm);
            add_contribution_to_emissions(*fbm, **di);

            likelihood_ += fbm->likelihood;
        }
        calculate_and_apply_new_parameters(hmm);

        LOG(INFO) << strprintf("%-3i  %-15.5g", i, likelihood_);
    }

    LOG(INFO) << hmm;

    return likelihood_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::setup(int num_states, int num_cols)
{
    if (fb_) delete fb_;
    fb_ = new ForwardBackwardAlgorithm<Alphabet_T, Subject_T>(*this);

    transitions_.clear();
    transitions_.resize(num_states + 1, num_states + 1);

    profiles_.clear();
    profiles_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>()));
    for (int k = 1; k <= num_states; ++k) {
        shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(num_cols));
        profiles_.push_back(profile_ptr);
    }

    likelihood_ = 0.0f;
    likelihood_prev_ = 0.0f;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_transitions(const ForwardBackwardMatrices& m,
                                                                               const HMM<Alphabet_T>& hmm)
{
    const int slen = m.f.num_rows() - 1;
    double prob_a_kl = 0.0f;
    for (const_transition_iterator ti = hmm.transitions_begin(); ti != hmm.transitions_end(); ++ti) {
        if (!transitions_.test(ti->from, ti->to)) transitions_[ti->from][ti->to] = 0.0f;

        if (ti->from == 0) {  // k=0 && l>0
            prob_a_kl = m.f[1][ti->to] * m.b[1][ti->to];
        } else {  // k>0 && l>0
            prob_a_kl = 0.0f;
            for (int i = 1; i < slen; ++i) {
                prob_a_kl += m.s[i+1] * m.f[i][ti->from] * m.b[i+1][ti->to] * ti->probability * m.e[i+1][ti->to];
            }
        }
        transitions_[ti->from][ti->to] = transitions_[ti->from][ti->to] + prob_a_kl;
        LOG(DEBUG1) << strprintf("tr[%i][%i] += %6.4f = %6.4f", ti->from, ti->to, prob_a_kl,
                                 static_cast<float>(transitions_[ti->from][ti->to]));
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                                                             const CountsProfile<Alphabet_T>& c)
{
    for (int k = 1; k < transitions_.num_rows(); ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];
        const int ci = p_k.center();

        for (int i = 1; i < m.f.num_rows(); ++i) {
            const int beg = std::max(0, i - ci - 1);
            const int end = std::min(c.num_cols() - 1, i + ci - 1);

            for(int h = beg; h <= end; ++h) {
                int j = h - i + ci + 1;
                for (int a = 0; a < p_k.alphabet_size(); ++a) {
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
    for (int k = 1; k < transitions_.num_rows(); ++k) {
        ContextProfile<Alphabet_T>& p_k = *profiles_[k];
        const int ci = p_k.center();

        for (int i = 1; i < m.f.num_rows(); ++i) {
            const int beg = std::max(0, i - ci - 1);
            const int end = std::min(s.length() - 1, i + ci - 1);

            for(int h = beg; h <= end; ++h) {
                int j = h - i + ci + 1;
                p_k[j][s[h]] += m.f[i][k] * m.b[i][k];
            }
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void BaumWelchTraining<Alphabet_T, Subject_T>::calculate_and_apply_new_parameters(HMM<Alphabet_T>& hmm)
{
    const int num_states = transitions_.num_rows() - 1;

    // Calculate and assign new emission probabilities
    for (int k = 1; k <= num_states; ++k) {
        normalize(*profiles_[k]);
        profiles_[k]->transform_to_logspace();
        hmm[k] = *profiles_[k];
        profiles_[k]->transform_to_linspace();
        reset(*profiles_[k], 0.0f);
    }

    // Calculate and assign new transition probabilities
    for (int k = 0; k <= num_states; ++k) {
        float sum = 0.0f;
        for (int l = 0; l <= num_states; ++l)
            if (transitions_.test(k,l)) sum += transitions_[k][l];

        if (sum != 0.0f) {
            float fac = 1.0f / sum;
            for (int l = 0; l <= num_states; ++l)
                if (transitions_.test(k,l)) {
                    hmm(k,l) = transitions_[k][l] * fac;
                    transitions_.erase(k,l);
                }
        } else if (k > 0) {
            // no out-transitions for this state -> connect to END-state
            transitions_[k][0] = 1.0f;
        } else {
            throw Exception("Unable calculate new HMM transition probabilities: BEGIN state has no out-transitions!");
        }
    }
}

}  // cs

#endif
