// Copyright 2009, Andreas Biegert

#ifndef SRC_TRANSITION_H_
#define SRC_TRANSITION_H_

namespace cs {

// Simple struct for transitions between HMM states.
struct Transition {
  // Simple Constructors
  Transition() : from(0), to(0), probability(0.0f) {}
  Transition(int f, int t, float p) : from(f), to(t), probability(p) {}
  ~Transition() {}

  operator float() const { return probability; }

  // Index of state from which the transition connects.
  const int from;
  // Index of state to which the transition connects.
  const int to;
  // Transition probability.
  float probability;
};  // Transition

struct AnchoredTransition {
  // Simple Constructors
  AnchoredTransition() : state(0), probability(0.0f) {}
  AnchoredTransition(int i, float p) : state(i), probability(p) {}
  ~AnchoredTransition() {}

  operator float() const { return probability; }

  // Index of state to/from which the transition connects.
  const int state;
  // Transition probability.
  float probability;
};  // AnchoredTransition

}  // namespace cs

#endif  // SRC_TRANSITION_H_
