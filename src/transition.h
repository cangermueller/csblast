// Copyright 2009, Andreas Biegert

#ifndef CS_TRANSITION_H_
#define CS_TRANSITION_H_

namespace cs {

// Simple struct for transitions between factor graph states.
struct Transition {
  // Default construction
  Transition() : source(0), target(0), prob(0.0) {}

  // Value construction
  Transition(int k, int l, double p) : source(k), target(l), prob(p) {}

  operator double() const { return prob; }

  // Index of source state of the transition
  unsigned int source;
  // Index of target state of the transition
  unsigned int target;
  // Transition weight.
  double prob;
};  // class Transition

// // Simple struct for transitions that are anchored at a state.
// struct AnchoredTransition {
//   AnchoredTransition() : state(0), weight(0.0f) {}
//   AnchoredTransition(int i, float p) : state(i), weight(p) {}
//   ~AnchoredTransition() {}

//   operator float() const { return weight; }

//   // Index of the source/target state of the transition.
//   const int state;
//   // Transition weight.
//   float weight;
// };  // class AnchoredTransition

}  // namespace cs

#endif  // CS_TRANSITION_H_
