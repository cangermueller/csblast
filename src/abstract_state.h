// Copyright 2009, Andreas Biegert

#ifndef SRC_ABSTRACT_STATE_H_
#define SRC_ABSTRACT_STATE_H_

namespace cs {

// A container class for profiles derived from alignments.
class AbstractState {
 public:
  virtual ~AbstractState() {}

  // Returns the index of the abstract state which acts as a unique identifier.
  virtual int index(int i) const = 0;
  // Sets the index of the abstract state.
  virtual int set_index(int i) const = 0;
};  // AbstractState

}  // namespace cs

#endif  // SRC_ABSTRACT_STATE_H_
