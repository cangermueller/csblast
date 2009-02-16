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

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class BaumWelchParams : public ForwardBackwardParams
{
  public:
    BaumWelchParams()
            : ForwardBackwardParams(),
              min_iterations_(3),
              max_iterations_(100),
              delta_likelihood_threshold_(1e-4)
    { }

    BaumWelchParams(const BaumWelchParams& params)
            : ForwardBackwardParams(params),
              min_iterations_(params.min_iterations_),
              max_iterations_(params.max_iterations_),
              delta_likelihood_(params.delta_likelihood_threshold_)
    { }

    float min_iterations() const { return min_iterations_; }
    float max_iterations() const { return max_iterations_; }
    float delta_likelihood_threshold() const { return delta_likelihood_threshold_; }
    void min_iterations(float min_iter) { return min_iterations_ = min_iter; }
    void max_iterations(float max_iter) { return max_iterations_ = max_iter; }
    void delta_likelihood(float threshold) { return delta_likelihood_threshold_ = threshold; }

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
    typedef typename std::vector< shared_ptr<Subject_T> > data_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;
    typedef shared_ptr<ForwardBackwardMatrices> fwd_bwd_matrices;

    BaumWelchTraining()
    { }

    BaumWelchTraining(const BaumWelchParams& params)
            : BaumWelchParams(params)
    { }

    float run(const HMM<Alphabet_T>& hmm, const data_vector& data)
    {
        LOG(DEBUG) << "Running Baum-Welch training on ...";
        LOG(DEBUG) << hmm;
        setup(hmm.num_states(), hmm[1].num_cols());

        for (int i = 1; i <= min_iterations() || i <= max_iterations() && delta_likelihood() > delta_likelihood_threshold(); ++i) {
            for (data_vector::const_iterator di = data.begin(); di != data.end(); ++di) {
                fwd_bwd_matrices fbm = fb->run(hmm, **di);
                add_contribution_to_transitions(*fbm, hmm);
                add_contribution_to_emissions(*fbm, **di);
            }
            calculate_new_parameters();
            apply_new_parameters();
        }

        return likelihood_;
    }

  private:
    void setup(int num_states, int num_cols)
    {
        fb = ForwardBackwardAlgorithm<Alphabet_T, Subject_T> fb(*this);

        transitions.clear();
        transitions_.resize(num_states + 1, num_states + 1);

        profiles_.clear();
        profiles_.push_back(shared_ptr());  // BEGIN-state has no emissions
        for (int k = 1; k <= num_states; ++k) {
            shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(num_cols));
            profiles_.push_back(profile_ptr);
        }

        likelihood_ = 0.0f;
        likelihood_prev = 0.0f;
    }

    void add_contribution_to_transitions(const fwd_bwd_matrices& m, const HMM& hmm)
    {
        for (int k = 0; k < transitions_.num_rows(); ++k) {
            for (int l = 0; l < transitions_.num_cols(); ++l) {
                if (k == 0 && l == 0 || !hmm.test_transition(k,l)) continue;
                if (!test(k,l)) transitions_[k][l] = 0.0f;

                for (int i = 1; i < m.f.num_rows()-1; ++i) {
                    float a_kl = m.c[i] * m.f[i][k] * m.b[i+1][l] * hmm(k,l) * m.e[i+1][l];
                    transitions_[k][l] = transitions_[k][l] + a_kl / m.likelihood;
                }
            }
        }
    }

    void add_contribution_to_emissions(const fwd_bwd_matrices& m, const CountsProfile& c)
    {
        for (int k = 1; k < transitions_.num_rows(); ++k) {
            ContextProfile& p_k = *profiles_[k];
            for (int j = 0; j < p_k.num_columns(); ++j) {
                for (int i = 1; i < m.f.num_rows(); ++i) {
                    for (int a = 0; a < p_k.alphabet_size(); ++a) {
                         p_k[j][a] += c[j][a] * m.f[i][k] * m.b[i][k] / m.likelihood;
                    }
                }
            }
        }
    }

    float delta_likelihood()
    {
        return likelihood_prev == 0.0f ? 1.0 : fabs(likelihood_ - likelihood_prev) / likelihood_prev;
    }

    // Instance of forward-backward algorithm to compute exptected number of transitions and emissions.
    ForwardBackwardAlgorithm* fb;
    // Transition counts
    sparse_matrix<float> transitions_;
    // Profiles with emission counts
    profiles_vector profiles_;
    // Likelihood of current iteration
    float likelihood_;
    // Likelihood of previous iteration
    float likelihood_prev_;
};

}  // cs

#endif
