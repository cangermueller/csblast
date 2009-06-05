// Copyright 2009, Andreas Biegert

#ifndef SRC_TRANSITION_H_
#define SRC_TRANSITION_H_

namespace cs {

// Simple struct for transitions between factor graph states.
struct Transition {
  Transition() : source(0), target(0), weight(0.0f) {}
  Transition(int f, int t, float p) : source(f), target(t), weight(p) {}
  ~Transition() {}

  operator float() const { return weight; }

  // Index of source state of the transition
  const int source;
  // Index of target state of the transition
  const int target;
  // Transition weight.
  float weight;
};  // class Transition

// Simple struct for transitions that are anchored at a state.
struct AnchoredTransition {
  AnchoredTransition() : state(0), weight(0.0f) {}
  AnchoredTransition(int i, float p) : state(i), weight(p) {}
  ~AnchoredTransition() {}

  operator float() const { return weight; }

  // Index of the source/target state of the transition.
  const int state;
  // Transition weight.
  float weight;
};  // class AnchoredTransition

}  // namespace cs

#endif  // SRC_TRANSITION_H_
