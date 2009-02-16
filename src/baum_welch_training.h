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
              max_iterations_(100),
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

    float run(const HMM<Alphabet_T>& hmm, const data_vector& data)
    {
        LOG(DEBUG) << "Running Baum-Welch training on ...";
        LOG(DEBUG) << hmm;
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
            //calculate_new_parameters();
            //apply_new_parameters();
            break;
        }

        for (int k = 1; k < hmm.num_states(); ++k) {
            normalize(*profiles_[k]);
            LOG(DEBUG) << *profiles_[k];
        }

        for (int k = 0; k <= hmm.num_states(); ++k) {
            for (int l = 0; l <= hmm.num_states(); ++l)
                std::cerr << strprintf("%6.4f  ", (float)transitions_[k][l]);
            std::cerr << std::endl;
        }

        return likelihood_;
    }

  private:
    void setup(int num_states, int num_cols)
    {
        if (fb_) delete fb_;
        fb_ = new ForwardBackwardAlgorithm<Alphabet_T, Subject_T>(*this);

        transitions_.clear();
        transitions_.resize(num_states + 1, num_states + 1);

        profiles_.clear();
        profiles_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>()));  // BEGIN-state has no emissions
        for (int k = 1; k <= num_states; ++k) {
            shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(num_cols));
            profiles_.push_back(profile_ptr);
        }

        likelihood_ = 0.0f;
        likelihood_prev_ = 0.0f;
    }

    void add_contribution_to_transitions(const ForwardBackwardMatrices& m, const HMM<Alphabet_T>& hmm)
    {
        for (int k = 0; k < transitions_.num_rows(); ++k) {
            for (int l = 0; l < transitions_.num_cols(); ++l) {
                if (k == 0 && l == 0 || !hmm.test_transition(k,l)) continue;
                if (!transitions_.test(k,l)) transitions_[k][l] = 0.0f;

                for (int i = 1; i < m.f.num_rows()-1; ++i) {
                    float a_kl = m.c[i] * m.f[i][k] * m.b[i+1][l] * hmm(k,l) * m.e[i+1][l];
                    transitions_[k][l] = transitions_[k][l] + a_kl / m.likelihood;
                }
            }
        }
    }

    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const CountsProfile<Alphabet_T>& c)
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
                        p_k[j][a] += c[h][a] * m.f[i][k] * m.b[i][k] / m.likelihood;
                    }
                }
            }
        }
    }

    void add_contribution_to_emissions(const ForwardBackwardMatrices& m, const Sequence<Alphabet_T>& s)
    {
        for (int k = 1; k < transitions_.num_rows(); ++k) {
            ContextProfile<Alphabet_T>& p_k = *profiles_[k];
            const int ci = p_k.center();

            for (int i = 1; i < m.f.num_rows(); ++i) {
                const int beg = std::max(0, i - ci - 1);
                const int end = std::min(s.length() - 1, i + ci - 1);

                for(int h = beg; h <= end; ++h) {
                    int j = h - i + ci + 1;
                    p_k[j][s[h]] += m.f[i][k] * m.b[i][k] / m.likelihood;
                }
            }
        }
    }

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
    float likelihood_;
    // Likelihood of previous iteration
    float likelihood_prev_;
};

}  // cs

#endif
