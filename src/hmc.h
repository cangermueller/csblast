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

#ifndef CS_HMC_H_
#define CS_HMC_H_

#include "crf-inl.h"
#include "func.h"
#include "parallel_tempering.h"
#include "progress_bar.h"
#include "sgd.h"

namespace cs {

// Parameter wrapper for various Hmc parameters
struct HmcParams {
  HmcParams()
      : nblocks(1),
        lfsteps(100),
        epsilon(1e-3),
        temp(1.0),
        epsilon_up(1.05),
        epsilon_down(0.8),
        alpha(1.3),
        sgd_epochs(20),
        seed(0) {}

  size_t nblocks;        // number of training blocks
  size_t lfsteps;        // number of leapfrog steps
  double epsilon;        // time step in leapfrog algorithm.
  double temp;           // base temperature for optimization
  double epsilon_up;     // epsilon up-scale after Metropolis acceptance
  double epsilon_down;   // epsilon down-scale after Metropolis rejection
  double alpha;          // geometric increase of temperatures in adjacent replicas
  int  sgd_epochs;       // number of SGD epochs in basin hopping
  unsigned int seed;     // seed for random number generators
};


template<class Abc>
struct HmcState : public DerivCrfFuncIO<Abc> {
  HmcState(const Crf<Abc>& c) : DerivCrfFuncIO<Abc>(c), steps(0) {}

  size_t steps;  // number of leapfrog steps already performed
};


template<class Abc, class TrainingPair>
struct LeapfrogProposal {
  LeapfrogProposal(const DerivCrfFunc<Abc, TrainingPair>& f,
                   const HmcParams& params,
                   int rank = 0)
      : func(f),
        nsteps(params.lfsteps),
        epsilon(params.epsilon),
        temp(params.temp * pow(params.alpha, rank)),
        nblocks(params.nblocks),
        epsilon_up(params.epsilon_up),
        epsilon_down(params.epsilon_down),
        gauss(0, 1, params.seed + rank) {}

  virtual ~LeapfrogProposal() {}

  // Given state 's1', sets state 's2' to a proposed candidate by leapfrog
  // integrating for 'nteps' time steps over energy landscape and returning the
  // acceptance probability of the new solution and rescale epsilon accordingly.
  virtual double operator() (HmcState<Abc>& s1,
                             HmcState<Abc>& s2,
                             ProgressBar* prog_bar = NULL) {
    assert_eq(s1.crf.nweights(), s2.crf.nweights());
    assert_eq((size_t)0, nsteps % (2 * nblocks));

    // Initialize gradient if this is the first leapfrog call
    if (s1.steps == 0) {
      if (prog_bar) prog_bar->Init((nsteps / nblocks + 1) *
                                   func.trainset.size() * s1.crf.size());
      InitGradient(s1, prog_bar);
    } else if (prog_bar) {
      prog_bar->Init((nsteps / nblocks) * func.trainset.size() * s1.crf.size());
    }
    s2 = s1;
    // Run leapfrog integration
    double ham_diff = LeapfrogIntegration(s2, prog_bar);
    // Rescale epsilon based on Metropolis acceptance probability
    double alph = MIN(1.0, exp(ham_diff));
    epsilon *= alph * epsilon_up + (1.0 - alph) * epsilon_down;

    return alph;
  }

  // Calculates initial gradient if this is the first leapfrog step
  void InitGradient(HmcState<Abc>& s, ProgressBar* prog_bar = NULL) const {
    // Initialize gradient, likelihood and prior
    s.loglike = 0.0;
    s.prior   = 0.0;
    Assign(s.grad_loglike, 0.0);
    Assign(s.grad_prior, 0.0);
    // Step backwards through training blocks to ensure that in the end s contains
    // gradient of first block, which is the one we will start with in leapfrog.
    for (int b = nblocks - 1; b >= 0; --b ) func.df(s, b, nblocks, prog_bar);
  }

  // Runs 'nsteps' leapfrog integration starting at state 's' and changing it in
  // place. Return value is the Metropolis acceptance probability.
  double LeapfrogIntegration(HmcState<Abc>& s, ProgressBar* prog_bar = NULL) {
    // Sample momentum terms from isotropic gaussian with stddev 1
    Vector<double> mom(s.crf.nweights());
    for (size_t i = 0; i < mom.size(); ++i) mom[i] = gauss();

    // Hamiltonian before leapfrog integration
    double ham1 = GetPotentialEnergy(s.loglike, s.prior) + GetKineticEnergy(mom);

    // First half step
    for (size_t i = 0; i < mom.size(); ++i) {
      if (i == 121) {
        LOG(DEBUG1) << strprintf("m[%4d]=%10.5f  ll[%4d]=%10.5f  pr[%4d]=%12.5f",
                                i, mom[i], i, s.grad_loglike[i] / temp, i,
                                s.grad_prior[i]);
      }
      mom[i] += 0.5 * epsilon * ((s.grad_loglike[i] / temp) + s.grad_prior[i]);
    }

    // Now we carry out all but the last leapfrog steps as full-steps
    size_t b = 0, c = 0;
    for (size_t step = 0; step < nsteps; ++step) {
      // Full-step update of CRF weights based on momenta
      for (size_t k = 0, i = 0; k < s.crf.size(); ++k) {
        s.crf[k].bias_weight += epsilon * mom[i++];
        for (size_t j = 0; j < s.crf.wlen(); ++j)
          for (size_t a = 0; a < Abc::kSize; ++a)
            s.crf[k].context_weights[j][a] += epsilon * mom[i++];
        for (size_t a = 0; a < Abc::kSize; ++a)
          s.crf[k].pc_weights[a] += epsilon * mom[i++];
        // Update pseudocounts based on new pseudocount weights
        UpdatePseudocounts(s.crf[k]);
      }
      // Reset likelihood at start of new scan through training blocks
      if (step % nblocks == 0) {
        s.loglike = 0.0;
        s.prior = 0.0;
      }
      // Calculate new gradient
      c = step % (2 * nblocks);  // index in forward-backward scan
      b = c < nblocks ? c : 2 * nblocks - c - 1;
      func.df(s, b, nblocks, prog_bar);

      // Update momenta with new gradient (half-step in last leapfrog step)
      const double eps = epsilon * (step < nsteps - 1 ? 1.0 : 0.5);
      for (size_t i = 0; i < mom.size(); ++i) {
        if (i == 121) {
          LOG(DEBUG1) << strprintf("m[%4d]=%10.5f  ll[%4d]=%10.5f  pr[%4d]=%12.5f",
                                  i, mom[i], i, s.grad_loglike[i] / temp, i,
                                  s.grad_prior[i]);
        }
        mom[i] += eps * ((s.grad_loglike[i] / temp) + s.grad_prior[i]);
      }
      s.steps++;
    }

    // Hamiltonian after leapfrog integration
    double ham2 = GetPotentialEnergy(s.loglike, s.prior) + GetKineticEnergy(mom);

    return ham1 - ham2;
  }

  double GetPotentialEnergy(double loglike, double prior) const {
    return (-loglike / temp) - prior;
  }

  double GetKineticEnergy(const Vector<double>& mom) const {
    double rv = 0.0;
    for (size_t i = 0; i < mom.size(); ++i) rv += SQR(mom[i]);
    return 0.5 * rv;
  }

  DerivCrfFunc<Abc, TrainingPair> func;  // functor for gradient calculation
  size_t nsteps;           // number of leapfrog steps
  double epsilon;          // time step epsilon
  double temp;             // temperature level
  size_t nblocks;          // number of training blocks
  double epsilon_up;       // epsilon up-scale after Metropolis acceptance
  double epsilon_down;     // epsilon down-scale after Metropolis rejection
  mutable Gaussian gauss;  // needed for sampling of momenta
};


template<class Abc, class TrainingPair>
struct BasinHoppingProposal : public LeapfrogProposal<Abc, TrainingPair> {

  BasinHoppingProposal(const DerivCrfFunc<Abc, TrainingPair>& f,
                       const HmcParams& hmc_params,
                       const SgdParams& sgd_params,
                       int rank = 0)
      : LeapfrogProposal<Abc, TrainingPair>(f, hmc_params, rank),
        sgd(f, sgd_params),
        sgd_epochs(hmc_params.sgd_epochs) {}

  virtual ~BasinHoppingProposal() {}

  // Given state 's1', sets state 's2' to a proposed candidate by leapfrog
  // integrating for 'nteps' time steps over energy landscape and returning the
  // acceptance probability of the new solution and rescale epsilon accordingly.
  virtual double operator() (HmcState<Abc>& s1,
                             HmcState<Abc>& s2,
                             ProgressBar* prog_bar = NULL) {
    assert_eq(s1.crf.nweights(), s2.crf.nweights());
    assert_eq((size_t)0, nsteps % (2 * nblocks));
    assert_ne(0, sgd_epochs);

    // Save potential energy before basin hopping
    double epot1 = GetPotentialEnergy(s1.loglike, s1.prior);
    // Initialize gradient if this is the first leapfrog call
    if (s1.steps == 0) {
      prog_bar->Init((nsteps / nblocks + sgd_epochs + 2) *
                     func.trainset.size() * s1.crf.size());
      InitGradient(s1, prog_bar);
    } else if (prog_bar) {
      prog_bar->Init((nsteps / nblocks + sgd_epochs + 1) *
                     func.trainset.size() * s1.crf.size());
    }
    s2 = s1;
    // Run leapfrog integration
    double ham_diff = LeapfrogIntegration(s2, prog_bar);
    // Rescale epsilon based on Metropolis acceptance probability
    double alph = MIN(1.0, exp(ham_diff));
    epsilon *= alph * epsilon_up + (1.0 - alph) * epsilon_down;
    // Save leapfrog proposal as new start point no matter what
    s1 = s2;
    LOG(ERROR) << "HMC loglike = " << s2.loglike << " HMC prior = " << s2.prior;
    // Run SGD to move proposed CRF into next local optimum
    SgdState<Abc> s(s2.crf);
    for (int epoch = 0; epoch < sgd_epochs; ++epoch) {
      sgd(s, prog_bar);
      LOG(ERROR) << "SGD loglike = " << s.loglike << " SGD prior = " << s.prior;
    }
    // Propose local optimum as new solution and calc its gradient and likelihhod
    s2.crf = s.crf;
    InitGradient(s2, prog_bar);
    // Calculate new potential energy
    double epot2 = GetPotentialEnergy(s2.loglike, s2.prior);

    LOG(ERROR) << "epot_diff = " << epot1 - epot2;

    return MIN(1.0, exp(epot1 - epot2));;
  }

  using LeapfrogProposal<Abc, TrainingPair>::func;
  using LeapfrogProposal<Abc, TrainingPair>::nsteps;
  using LeapfrogProposal<Abc, TrainingPair>::nblocks;
  using LeapfrogProposal<Abc, TrainingPair>::epsilon;
  using LeapfrogProposal<Abc, TrainingPair>::epsilon_up;
  using LeapfrogProposal<Abc, TrainingPair>::epsilon_down;
  Sgd<Abc, TrainingPair> sgd;  // SGD algorithm encapsulation
  int sgd_epochs;              // number of SGD epochs to run after leapfrog
};


// Driver function for performin 'm' steps of HMC sampling, starting at state 's'
// and proposing new candaidate solutions with proposal functor 'propose'
template<class Abc, class TrainingPair>
double HmcStep(size_t m,                                         // sampling steps
               HmcState<Abc>& s,                                 // state object
               LeapfrogProposal<Abc, TrainingPair>& propose,     // proposal functor
               const CrfFunc<Abc, TrainingPair>& loglike, // LL on valid. set
               unsigned int seed = 0,
               ParallelTempering<Abc, TrainingPair>* parallel_tempering = NULL,
               FILE* fout = stdout) {
  HmcState<Abc> sprop(s);     // storage for next sample
  HmcState<Abc> sbest(s);     // best solution found so far
  int accept = 0;             // number of accepted samples
  bool is_accepted = false;   // was last sample accepted?
  Ran ran(seed);
  scoped_ptr<ProgressBar> prog_bar;
  double alph, ll, ll_best;

  if (fout) {
    prog_bar.reset(new ProgressBar(fout, 12));
    fprintf(fout, "%-4s %5s %7s  %-12s %8s  %9s  %5s  %3s  %5s  %7s\n", "No",
            "Temp", "Epsilon", "Proposal", "LL-Train", "Prior", "Prob", "y/n",
            "Acc", "LL-Val");
    fprintf(fout, "%s\n", std::string(80, '-').c_str());
  }

  ll_best = loglike(sprop.crf);
  for (size_t i = 0; i < m; ++i) {
    if (fout) {
      fprintf(fout, "%-4zu %5.2f  %6.1g  ", i+1, propose.temp, propose.epsilon);
      fflush(fout);
    }

    // Calculate new proposal and its likelihood on the validation set
    alph = propose(s, sprop, prog_bar.get());
    ll = loglike(sprop.crf);
    // Evaluate metripolis criterion
    is_accepted = (ran.doub() < alph);
    // Update state based on Metropolis decision
    if (is_accepted) {
      s = sprop;
      accept++;
    }
    // Keep track of best solution on validation set
    if (ll > ll_best) {
      ll_best = ll;
      sbest = sprop;
    }

    if (fout) {
      fprintf(fout, " %8.4f %10.4f  %5.1f  %3s  %5.1f %8.4f\n",
              sprop.loglike / propose.func.trainset.size(),
              sprop.prior / sprop.crf.nweights(), alph * 100,
              is_accepted ? "yes" : "no", 100 * accept / (i + 1.0),
              ll / loglike.trainset.size());
    }

    if (parallel_tempering)
      parallel_tempering->ReplicaExchange(s.loglike, propose.temp);
  }
  s = sbest;
  return accept / static_cast<double>(m);
}

}  // namespace cs

#endif  // CS_HMC_H_
