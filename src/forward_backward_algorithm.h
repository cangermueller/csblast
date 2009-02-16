#ifndef CS_FORWARD_BACKWARD_ALGORITHM_H
#define CS_FORWARD_BACKWARD_ALGORITHM_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Forward-Backward algorithm encapsulation.

#include <valarray>

#include "context_profile.h"
#include "hmm.h"
#include "log.h"
#include "profile_matcher.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

class ForwardBackwardParams
{
  public:
    ForwardBackwardParams()
            : ignore_begin_transitions_(false),
              ignore_end_transitions_(false),
              ignore_profile_context_(false),
              weight_center_(1.6f),
              weight_decay_(0.85f)
    { }

    ForwardBackwardParams(const ForwardBackwardParams& p)
            : ignore_begin_transitions_(p.ignore_begin_transitions_),
              ignore_end_transitions_(p.ignore_end_transitions_),
              ignore_profile_context_(p.ignore_profile_context_),
              weight_center_(p.weight_center_),
              weight_decay_(p.weight_decay_)
    { }

    ForwardBackwardParams& ignore_begin_transitions(bool f) { ignore_begin_transitions_ = f; return *this; }
    ForwardBackwardParams& ignore_end_transitions(bool f) { ignore_end_transitions_ = f; return *this; }
    ForwardBackwardParams& ignore_profile_context(bool f) { ignore_profile_context_ = f; return *this; }
    ForwardBackwardParams& weight_center(float w) { weight_center_ = w; return *this; }
    ForwardBackwardParams& weight_decay(float b) { weight_decay_ = b; return *this; }

    bool ignore_begin_transitions() const { return ignore_begin_transitions_; }
    bool ignore_end_transitions() const { return ignore_end_transitions_; }
    bool ignore_profile_context() const { return ignore_profile_context_; }
    float weight_center() const { return weight_center_; }
    float weight_decay() const { return weight_decay_; }

  private:
    bool ignore_begin_transitions_;
    bool ignore_end_transitions_;
    bool ignore_profile_context_;
    float weight_center_;
    float weight_decay_;
};

struct ForwardBackwardMatrices
{
    typedef float value_type;

    ForwardBackwardMatrices()
            : likelihood(0.0)
    { }

    ForwardBackwardMatrices(int nrows, int ncols)
            : f(nrows, ncols, 0.0),
              b(nrows, ncols, 0.0),
              e(nrows, ncols, 1.0),
              c(0.0, nrows),
              likelihood(0.0)
    { }

    friend std::ostream& operator<< (std::ostream& out, const ForwardBackwardMatrices& m)
    {
        out << "Forward matrix f[i][k]:" << std::endl;
        for (int i = 1; i < m.f.num_rows(); ++i) {
            for (int k = 1; k < m.f.num_cols(); ++k) {
                out << strprintf("%6.4f  ", m.f[i][k]);
            }
            out << std::endl;
        }
        out << "Backward matrix b[i][k]:" << std::endl;
        for (int i = 1; i < m.f.num_rows(); ++i) {
            for (int k = 1; k < m.f.num_cols(); ++k) {
                out << strprintf("%6.4f  ", m.b[i][k]);
            }
            out << std::endl;
        }
        out << "Forward row-sums before scaling c[i]:" << std::endl;
        for (int i = 1; i < m.f.num_rows(); ++i) {
            out << strprintf("%6.2f  ", m.c[i]);
        }
        out << "\nEmission probabilities e[i][k]:" << std::endl;
        for (int i = 1; i < m.f.num_rows(); ++i) {
            for (int k = 1; k < m.f.num_cols(); ++k) {
                out << strprintf("%6.4f  ", m.e[i][k]);
            }
            out << std::endl;
        }
        out << "Posterior probability matrix p[i][k]:" << std::endl;
        for (int i = 1; i < m.f.num_rows(); ++i) {
            for (int k = 1; k < m.f.num_cols(); ++k) {
                value_type p = 100.0 * m.f[i][k] * m.b[i][k] / m.likelihood;
                out << strprintf("%6.2f  ", p);
            }
            out << std::endl;
        }
        out << strprintf("Likelihood = %-7.2g\n", m.likelihood);
        return out;
    }

    // forward matrix
    matrix<value_type> f;
    // backward matrix
    matrix<value_type> b;
    // emission probability matrix
    matrix<value_type> e;
    // forward matrix column sum c(n) before normalization
    std::valarray<value_type> c;
    // likelihood P(x)
    value_type likelihood;
};

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardAlgorithm : public ForwardBackwardParams
{
  public:
    typedef typename ForwardBackwardMatrices::value_type value_type;

    ForwardBackwardAlgorithm()
    { }

    ForwardBackwardAlgorithm(const ForwardBackwardParams& params)
            : ForwardBackwardParams(params)
    { }

    shared_ptr<ForwardBackwardMatrices> run(const HMM<Alphabet_T>& hmm, const Subject_T<Alphabet_T>& subject)
    {
        LOG(DEBUG) << "Running forward-backward algorithm on ...";
        LOG(DEBUG) << hmm;
        LOG(DEBUG) << subject;

        ProfileMatcher<Alphabet_T> matcher;
        if (!ignore_profile_context() && hmm[1].num_cols() > 1)
            matcher.init_weights(hmm[1].num_cols(), weight_center(), weight_decay());

        shared_ptr<ForwardBackwardMatrices> matrices( new ForwardBackwardMatrices(subject.length() + 1,
                                                                                  hmm.num_states() + 1) );

        matrices->likelihood = forward(hmm, subject, matcher, *matrices);
        backward(hmm, subject, *matrices);

        return matrices;
    }

  private:
    typedef typename HMM<Alphabet_T>::const_state_iterator const_state_iterator;
    typedef typename State<Alphabet_T>::const_transition_iterator const_transition_iterator;

    value_type forward(const HMM<Alphabet_T>& hmm,
                       const Subject_T<Alphabet_T>& subject,
                       const ProfileMatcher<Alphabet_T>& matcher,
                       ForwardBackwardMatrices& m)
    {
        // TODO: pay attention to ignore BEGIN/END state flags.
        LOG(DEBUG) << "Forward pass ...";

        const int length = subject.length();
        // Initialization
        m.f[0][0] = 1.0;
        for (int i = 1; i <= length; ++i) {
            LOG(DEBUG) << strprintf("i=%i", i);
            m.c[i] = 0.0;

            for (const_state_iterator s_l = hmm.states_begin(); s_l != hmm.states_end(); ++s_l) {
                value_type f_il = 0.0;
                LOG(DEBUG1) << strprintf("f[%i][%i] = 0", i, (**s_l).index());

                for (const_transition_iterator t_kl = (**s_l).in_transitions_begin(); t_kl != (**s_l).in_transitions_end(); ++t_kl) {
                    f_il += m.f[i-1][t_kl->state] * t_kl->probability;

                    LOG(DEBUG2) << strprintf("f[%i][%i] += f[%i][%i]=%-7.2e * tr[%i][%i]=%-7.5f",
                                             i, (**s_l).index(), i-1, t_kl->state, m.f[i-1][t_kl->state],
                                             t_kl->state, (**s_l).index(), t_kl->probability);
                }

                m.e[i][(**s_l).index()] = matcher(**s_l, subject, i-1);
                f_il *= m.e[i][(**s_l).index()];
                LOG(DEBUG2) << strprintf("f[%i][%i] *= e[%i][%i]=%-7.5f", i, (**s_l).index(), i, (**s_l).index(), m.e[i][(**s_l).index()]);

                m.f[i][(**s_l).index()] = f_il;
                m.c[i] += f_il;
                LOG(DEBUG1) << strprintf("f[%i][%i] = %-7.2e", i, (**s_l).index(), f_il);
            }

            // Rescale forward values
            value_type scale_fac = 1.0 / m.c[i];
            for (int l = 1; l <= hmm.num_states(); ++l)
                m.f[i][l] *= scale_fac;
        }

        value_type likelihood = 0.0;
        for (const_transition_iterator t_k0 = hmm[0].in_transitions_begin(); t_k0 != hmm[0].in_transitions_end(); ++t_k0)
            likelihood += m.f[length][t_k0->state] * t_k0->probability;

        LOG(DEBUG) << strprintf("Forward likelihood = %-7.2g", likelihood);
        return likelihood;
    }

    value_type backward(const HMM<Alphabet_T>& hmm,
                        const Subject_T<Alphabet_T>& subject,
                        ForwardBackwardMatrices& m)
    {
        // TODO: pay attention to ignore BEGIN/END state flags.
        LOG(DEBUG) << "Backward pass ...";

        const int length = subject.length();
        // Initialization
        for (const_transition_iterator t_k0 = hmm[0].in_transitions_begin(); t_k0 != hmm[0].in_transitions_end(); ++t_k0)
            m.b[length][t_k0->state] = t_k0->probability;

        for (int i = length-1; i >= 1; --i) {
            LOG(DEBUG) << strprintf("i=%i", i);
            for (const_state_iterator s_k = hmm.states_begin(); s_k != hmm.states_end(); ++s_k) {
                value_type b_ik = 0.0;
                LOG(DEBUG1) << strprintf("b[%i][%i] = 0", i, (**s_k).index());

                for (const_transition_iterator t_kl = (**s_k).out_transitions_begin(); t_kl != (**s_k).out_transitions_end(); ++t_kl) {
                    if (t_kl->state == 0) continue;  // ignore transitions to END-state
                    b_ik += t_kl->probability * m.e[i+1][t_kl->state] * m.b[i+1][t_kl->state];

                    LOG(DEBUG2) << strprintf("b[%i][%i] += tr[%i][%i]=%-7.5f * e[%i][%i]=%-7.5f * f[%i][%i]=%-7.2e",
                                             i, (**s_k).index(), t_kl->state, (**s_k).index(), t_kl->probability,
                                             i+1, t_kl->state, m.e[i+1][t_kl->state], i+1, t_kl->state, m.b[i+1][t_kl->state]);
                }
                m.b[i][(**s_k).index()] = b_ik / m.c[i+1];
                LOG(DEBUG1) << strprintf("b[%i][%i] = %-7.2e", i, (**s_k).index(), b_ik);
            }
        }

        value_type likelihood = 0.0;
        for (const_transition_iterator t_0l = hmm[0].out_transitions_begin(); t_0l != hmm[0].out_transitions_end(); ++t_0l)
            likelihood += t_0l->probability * m.e[1][t_0l->state] * m.b[1][t_0l->state] / m.c[1];

        LOG(DEBUG) << strprintf("Backward likelihood = %-7.2g", likelihood);
        return likelihood;
    }

};

}  // cs

#endif
