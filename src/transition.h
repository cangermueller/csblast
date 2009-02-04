#ifndef CS_STATE_H
#define CS_STATE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The state object of a context HMM.

//#include <iostream>

#include "context_profile.h"
#include "exception.h"
#include "context_profile.h"
#include "shared_ptr.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class State;

struct Transition
{
    // Simple Constructors
    Transition() : from(0), to(0), probability(0.0f) {}
    Transition(int f, int t, float p) : from(f), to(t), probability(p) {}
    ~Transition() {}

    // Index of state from which the transition connects.
    int from;
    // Index of state to which the transition connects.
    int to;
    // Transition probability.
    float probability;
};  // Transition



// Functor predicate that returns true if a transition has a transition probability below given threshold.
class FailsProbabilityThreshold : public std::unary_function<Transition, bool>
{
  public:
    FailsProbabilityThreshold(float threshold) : threshold_(threshold) {}

    bool operator() (const Transition& transition) const
    {
        return transition.probability < threshold_;
    }

  private:
    const float threshold_;
};

// Functor predicate that returns true if a transition transits from or to a given state.
class TransitsState : public std::unary_function<Transition, bool>
{
  public:
    TransitsState(int state) : state_(state) {}

    bool operator() (const Transition& transition) const
    {
        return transition.from == state_;
    }

  private:
    const int state_;
};

// // Functor predicate that returns true if a transition transits from or to a given state.
// class TransitsState : public std::unary_function<Transition, bool>
// {
//   public:
//     TransitsState(int state) : state_(state) {}

//     bool operator() (const Transition& transition) const
//     {
//         return transition.state == state_;
//     }

//   private:
//     const int state_;
// };

}  // cs

#endif
