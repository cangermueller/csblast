#ifndef CS_TRANSITION_H
#define CS_TRANSITION_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Simple struct for transitions between HMM states.

namespace cs
{

struct Transition
{
    // Simple Constructors
    Transition() : from(0), to(0), probability(0.0f) {}
    Transition(int f, int t, float p) : from(f), to(t), probability(p) {}
    ~Transition() {}

    // Index of state from which the transition connects.
    const int from;
    // Index of state to which the transition connects.
    const int to;
    // Transition probability.
    float probability;
};  // Transition

}  // cs

#endif