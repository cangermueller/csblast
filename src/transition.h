/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
