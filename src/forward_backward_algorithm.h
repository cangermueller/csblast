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

// Forward declaration
template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardAlgorithm;

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardParameters
{
  public:
    ForwardBackwardParameters()
            : ignore_begin_transitions_(false),
              ignore_end_transitions_(false),
              ignore_profile_context_(false),
              weight_center_(1.6f),
              weight_decay_(0.85f)
    { }

    // sets all the default values for each data member
    ForwardBackwardParameters& ignore_begin_transitions() { ignore_begin_transitions_ = true; return *this; }
    ForwardBackwardParameters& ignore_end_transitions() { ignore_end_transitions_ = true; return *this; }
    ForwardBackwardParameters& ignore_profile_context() { ignore_profile_context_ = true; return *this; }
    ForwardBackwardParameters& weight_center(float w) { weight_center_ = w; return *this; }
    ForwardBackwardParameters& weight_decay(float b) { weight_decay_ = b; return *this; }

    //  private:
    friend class ForwardBackwardAlgorithm<Alphabet_T, Subject_T>;

    bool ignore_begin_transitions_;
    bool ignore_end_transitions_;
    bool ignore_profile_context_;
    float weight_center_;
    float weight_decay_;
};

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ForwardBackwardAlgorithm
{
  public:
    ForwardBackwardAlgorithm(const ForwardBackwardParameters<Alphabet_T, Subject_T>& params)
            : ignore_begin_transitions_(params.ignore_begin_transitions_),
              ignore_end_transitions_(params.ignore_end_transitions_)
              // TODO: create profile matcher with weight_center, weight_decay, etc.
    { }

    bool ignore_begin_transitions_;
    bool ignore_end_transitions_;
};

}  // cs

#endif
