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

#ifndef CS_HMM_H_
#define CS_HMM_H_

#include <google/sparsetable>

#include "context_profile.h"
#include "context_library.h"
#include "emission.h"
#include "pseudocounts-inl.h"

namespace cs {

// Forward declarations
template<class Abc>
class Hmm;

// Strategy class for initializing a context HMM
template<class Abc>
class HmmInit {
 public:
  HmmInit() {}
  virtual ~HmmInit() {}
  virtual void operator() (Hmm<Abc>& hmm) const = 0;
};


// Simple struct for transitions between context profiles.
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


// A container class of K context profiles and transitions probabilities between
// them. Alltogether representing the most common sequence motifs in a training
// database of proteins/DNA sequences.
template<class Abc>
class Hmm : public ContextLibrary<Abc> {
 public:
  typedef ContextProfile<Abc>* ProfileIter;
  typedef const ContextProfile<Abc>* ConstProfileIter;

  // Constructs an empty HMM of given dimenions.
  Hmm(size_t size, size_t wlen);

  // Constructs an HMM from serialized data read from input stream.
  explicit Hmm(FILE* fin);

  // Constructs an HMM with a specific init-strategy encapsulated by an
  // initializer.
  Hmm(size_t size, size_t wlen, const HmmInit<Abc>& init);

  // Nothing to do here
  virtual ~Hmm() {}

  // Writes the HMM in serialization format to output stream.
  void Write(FILE* fout) const;

 private:
   // Initializes the HMM from serialized data read from stream.
  void Read(FILE* fin);

  size_t wlen_;                             // size of context window.
  Vector< ContextProfile<Abc> > profiles_;  // context profiles ordered by index.
};  // ContextLibrary


// Prints the library in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const ContextLibrary<Abc>& lib) {
  out << "ContextLibrary" << std::endl;
  out << "size:\t" << lib.size() << std::endl;
  out << "wlen:\t" << lib.wlen() << std::endl;
  for (size_t k = 0; k < lib.size(); ++k) out << lib[k];
  return out;
}

// Transforms probabilites in context profiles to log-space and sets 'is_log' flag.
template<class Abc>
void TransformToLog(ContextLibrary<Abc>& lib);

// Transforms probabilites in context profiles to lin-space and sets 'is_log' flag.
template<class Abc>
void TransformToLin(ContextLibrary<Abc>& lib);

// Calculates posterior probs for a context library and sequence window X_i
// centered at index 'i' and writes them to array 'pp'. Caller is responsible for
// making sure that 'pp' has sufficient length. Return value is log sum of all
// individual emission terms.
template<class Abc, class ContextInput>
double CalculatePosteriorProbs(const ContextLibrary<Abc>& lib,
                               const Emission<Abc>& emission,
                               const ContextInput& input,
                               size_t i,
                               double* pp);

// Strategy for initializing library by sampling from training set of count
// profiles, optionally adding pseudocounts.
template<class Abc>
class SamplingLibraryInit : public LibraryInit<Abc> {
 public:
  typedef std::vector< CountProfile<Abc> > TrainingSet;

  SamplingLibraryInit(const TrainingSet& trainset,
                      const Pseudocounts<Abc>& pc,
                      const Admix& admix,
                      unsigned int seed = 0)
      : trainset_(trainset),
        pc_(pc),
        admix_(admix),
        seed_(seed) {}

  virtual ~SamplingLibraryInit() {}

  virtual void operator() (ContextLibrary<Abc>& lib) const;

 private:
  const TrainingSet& trainset_;
  const Pseudocounts<Abc>& pc_;
  const Admix& admix_;
  const unsigned int seed_;
};  // SamplingLibraryInit


// Strategy that initializes profile probs by sammpling from gaussian distribution
// with mean at background frequencies.
template<class Abc>
class GaussianLibraryInit : public LibraryInit<Abc> {
 public:
  GaussianLibraryInit(double sigma,
                      const SubstitutionMatrix<Abc>& sm,
                      unsigned int seed = 0)
      : sigma_(sigma), sm_(sm), seed_(seed) {}

  virtual ~GaussianLibraryInit() {}

  virtual void operator() (ContextLibrary<Abc>& lib) const;

 protected:
  double sigma_;
  const SubstitutionMatrix<Abc>& sm_;
  unsigned int seed_;
};  // class GaussianLibraryInit

}  // namespace cs

#endif  // CS_CONTEXT_LIBRARY_H_
