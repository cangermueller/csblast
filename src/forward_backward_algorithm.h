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
    ForwardBackwardMatrices()
            : p_forward(0.0f)
    { }

    ForwardBackwardMatrices(int nrows, int ncols)
            : f(nrows, ncols, 0.0f),
              b(nrows, ncols, 0.0f),
              p_forward(0.0f)
    { }

    // forward matrix
    matrix<float> f;
    // backward matrix
    matrix<float> b;
    // likelihood P(x)
    float p;
};

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardAlgorithm : public ForwardBackwardParams
{
  public:
    ForwardBackwardAlgorithm(const ForwardBackwardParams& params)
            : ForwardBackwardParams(params)
    { }

    shared_ptr<ForwardBackwardMatrices> run(const HMM& hmm, const Subject_T& subject)
    {
        const ProfileMatcher matcher(ignore_profile_context() ? 1 : hmm[0].num_cols(),
                                     weight_center(),
                                     weight_decay());
        shared_ptr<ForwardBackwardMatrices> matrices(subject.length() + 1, hmm.num_states() + 1);

        float p_forward  = forward(hmm, subject, matcher, matrices->f);
        float p_backward = backward(hmm, subject, matcher, matrices->b);
        matrices->p = p_forward;

        return matrices;
    }

  private:
    typedef typename HMM<Alphabet_T>::const_state_iterator const_state_iterator;
    typedef typename State<Alphabet_T>::const_transition_iterator const_transition_iterator;

    float forward(const HMM<Alphabet_T>& hmm,
                  const Subject_T<Alphabet_T>& subject,
                  const ProfileMatcher<Alphabet_T>& matcher,
                  matrix<float>& f)
    {
        const int length = subject.length();

        // TODO: pay attention to ignore BEGIN/END state flags.

        f[0][0] = 1.0;  // we start from the BEGIN state
        for (int i = 1; i <= length; ++i) {
            for (const_state_iterator s_l = hmm.states_begin(); s_l != hmm.states_end(); ++s_l) {
                float f_il = matcher.match(*s_l, subject, i);
                for (const_transition_iterator t_kl = s_l.in_transitions_begin(); t_kl != s_l.in_transitions_end(); ++t_kl) {
                    f_il += f[i-1][t_kl->state] * t_kl->probability;
                }
            }
        }

        float p_forward = 0.0f;
        for (const_transition_iterator t_k0 = hmm[0].in_transitions_begin(); t_k0 != hmm[0].in_transitions_end(); ++t_k0)
            p_forward += f[length][t_k0->state] * t_k0->probability;

        return p_forward;
    }

    float backward(const HMM<Alphabet_T>& hmm,
                   const Subject_T<Alphabet_T>& subject,
                   const ProfileMatcher<Alphabet_T>& matcher,
                   matrix<float>& b)
    {
        // TODO
    }

};

}  // cs

#endif
