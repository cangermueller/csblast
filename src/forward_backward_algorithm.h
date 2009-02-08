#ifndef CS_FORWARD_BACKWARD_ALGORITHM_H
#define CS_FORWARD_BACKWARD_ALGORITHM_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Forward-Backward algorithm encapsulation.

#include "context_profile.h"
#include "hmm.h"

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

    bool ignore_begin_transitions() { return ignore_begin_transitions_; }
    bool ignore_end_transitions() { return ignore_end_transitions_; }
    bool ignore_profile_context() { return ignore_profile_context_; }
    float weight_center() { return weight_center_; }
    float weight_decay() { return weight_decay_; }

  private:
    bool ignore_begin_transitions_;
    bool ignore_end_transitions_;
    bool ignore_profile_context_;
    float weight_center_;
    float weight_decay_;
};

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardAlgorithm : public ForwardBackwardParams
{
  public:
    ForwardBackwardAlgorithm(const ForwardBackwardParams& params)
            : ForwardBackwardParams(params)
    { }
};

}  // cs

#endif
