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
    typedef double value_type;

    ForwardBackwardMatrices()
            : p_forward(0.0)
    { }

    ForwardBackwardMatrices(int nrows, int ncols)
            : f(nrows, ncols, 0.0),
              b(nrows, ncols, 0.0),
              p_forward(0.0)
    { }

    // forward matrix
    matrix<value_type> f;
    // backward matrix
    matrix<value_type> b;
    // likelihood P(x)
    value_type p_forward;
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
        LOG(DEBUG) << "Running forward-backward algorithm ...";
        LOG(DEBUG) << hmm;
        LOG(DEBUG) << subject;

        ProfileMatcher<Alphabet_T> matcher;
        if (!ignore_profile_context() && hmm[1].num_cols() > 1)
            matcher.init_weights(hmm[1].num_cols(), weight_center(), weight_decay());

        shared_ptr<ForwardBackwardMatrices> matrices( new ForwardBackwardMatrices(subject.length() + 1,
                                                                                  hmm.num_states() + 1) );

        value_type p_forward  = forward(hmm, subject, matcher, matrices->f);
        LOG(DEBUG) << "P_forward=" << p_forward;
        backward(hmm, subject, matcher, matrices->b);

        matrices->p_forward = p_forward;
        return matrices;
    }

  private:
    typedef typename HMM<Alphabet_T>::const_state_iterator const_state_iterator;
    typedef typename State<Alphabet_T>::const_transition_iterator const_transition_iterator;

    value_type forward(const HMM<Alphabet_T>& hmm,
                  const Subject_T<Alphabet_T>& subject,
                  const ProfileMatcher<Alphabet_T>& matcher,
                  matrix<value_type>& f)
    {
        // TODO: pay attention to ignore BEGIN/END state flags.
        LOG(DEBUG) << "Forward pass ...";

        const int length = subject.length();
        f[0][0] = 1.0;  // we start from the BEGIN state
        for (int i = 1; i <= length; ++i) {
            LOG(DEBUG) << "i=" << i;
            for (const_state_iterator s_l = hmm.states_begin(); s_l != hmm.states_end(); ++s_l) {
                LOG(DEBUG1) << "Calculating f[i=" << i << "][l=" << (**s_l).index() << "] ...";
                value_type f_il = matcher.match(**s_l, subject, i-1);
                LOG(DEBUG2) << "f[i=" << i << "][l=" << (**s_l).index() << "] = " << f_il;
                for (const_transition_iterator t_kl = (**s_l).in_transitions_begin(); t_kl != (**s_l).in_transitions_end(); ++t_kl) {
                    LOG(DEBUG3) << "Processing tr[k=" << t_kl->state << "][l=" << (**s_l).index() << "] ...";
                    f_il += f[i-1][t_kl->state] * t_kl->probability;
                    LOG(DEBUG4) << "f[i=" << i << "][l=" << (**s_l).index() << "] += f[" << i-1 << "][" << t_kl->state << "]="
                                << f[i-1][t_kl->state] << " * " << "tr[" << t_kl->state << "][" << (**s_l).index() << "]="
                                << t_kl->probability;
                }
                f[i][(**s_l).index()] = f_il;
                LOG(DEBUG2) << "f[i=" << i << "][k=" << (**s_l).index() << "] = " << f_il;;
            }
        }

        value_type p_forward = 0.0;
        for (const_transition_iterator t_k0 = hmm[0].in_transitions_begin(); t_k0 != hmm[0].in_transitions_end(); ++t_k0)
            p_forward += f[length][t_k0->state] * t_k0->probability;

        return p_forward;
    }

    value_type backward(const HMM<Alphabet_T>& hmm,
                   const Subject_T<Alphabet_T>& subject,
                   const ProfileMatcher<Alphabet_T>& matcher,
                   matrix<value_type>& b)
    {
        // TODO: pay attention to ignore BEGIN/END state flags.
        LOG(DEBUG) << "Backward pass ...";

        const int length = subject.length();
        // Initialization
        for (const_transition_iterator t_k0 = hmm[0].in_transitions_begin(); t_k0 != hmm[0].in_transitions_end(); ++t_k0)
            b[length][t_k0->state] = t_k0->probability;

        for (int i = length-1; i >= 1; --i) {
            LOG(DEBUG) << "i=" << i;
            for (const_state_iterator s_k = hmm.states_begin(); s_k != hmm.states_end(); ++s_k) {
                LOG(DEBUG1) << "Calculating b[i=" << i << "][k=" << (**s_k).index() << "] ...";

                value_type b_ik = 0.0f;
                LOG(DEBUG2) << "b[i=" << i << "][k=" << (**s_k).index() << "] = " << b_ik;
                for (const_transition_iterator t_kl = (**s_k).out_transitions_begin(); t_kl != (**s_k).out_transitions_end(); ++t_kl) {
                    if (t_kl->state == 0) continue;  // ignore transitions to END-state
                    LOG(DEBUG3) << "Processing tr[k=" << (**s_k).index() << "][l=" << t_kl->state << "] ...";

                    b_ik += t_kl->probability * matcher.match(hmm[t_kl->state], subject, i) * b[i+1][t_kl->state];

                    LOG(DEBUG4) << "b[i=" << i << "][k=" << (**s_k).index() << "] += tr[" << (**s_k).index()
                                << "][" << t_kl->state << "]=" << t_kl->probability << " * e[" << t_kl->state
                                << "][" << i << "]=" << matcher.match(hmm[t_kl->state], subject, i) << "* b[i="
                                << i+1 << "][" << t_kl->state << "]=" << b[i+1][t_kl->state];
                }
                b[i][(**s_k).index()] = b_ik;
                LOG(DEBUG2) << "b[i=" << i << "][k=" << (**s_k).index() << "] = " << b_ik;
            }
        }

        value_type p_backward = 0.0;
        for (const_transition_iterator t_0l = hmm[0].out_transitions_begin(); t_0l != hmm[0].out_transitions_end(); ++t_0l)
            p_backward += t_0l->probability * matcher.match(hmm[t_0l->state], subject, 0) * b[1][t_0l->state];

        return p_backward;
    }

};

}  // cs

#endif
