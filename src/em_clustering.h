// Copyright 2009, Andreas Biegert

#ifndef CS_EM_CLUSTERING_H_
#define CS_EM_CLUSTERING_H_

#include "context_library-inl.h"
#include "emission.h"
#include "progress_bar.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Parameter wrapper for various Hmc parameters
struct EMClusteringParams {
  EMClusteringParams()
      : weight_center(1.6),
        weight_decay(0.85),
        pca(0.01) {}

  double weight_center;
  double weight_decay;
  double pca;
};

template<class Abc>
struct EMClustering {
  typedef std::vector<CountProfile<Abc> > TrainingSet;

  EMClustering(const TrainingSet& tset,
               const ContextLibrary<Abc>& cl,
               const SubstitutionMatrix<Abc>& m,
               double w_center,
               double w_decay,
               double a)
      : trainset(tset),
        sm(m),
        lib(cl),
        priors(lib.size(), 0.0),
        probs(lib.size(), Profile<Abc>(lib.wlen(), 0.0)),
        weight_center(w_center),
        weight_decay(w_decay),
        pca(a),
        loglike(-DBL_MAX) {}

  EMClustering(const TrainingSet& tset,
               const ContextLibrary<Abc>& cl,
               const SubstitutionMatrix<Abc>& m,
               const EMClusteringParams& params)
      : trainset(tset),
        sm(m),
        lib(cl),
        priors(lib.size(), 0.0),
        probs(lib.size(), Profile<Abc>(lib.wlen(), 0.0)),
        weight_center(params.weight_center),
        weight_decay(params.weight_decay),
        pca(params.pca),
        loglike(-DBL_MAX) {}

  double EStep(ProgressBar* prog_bar = NULL) {
    const int ntrain = trainset.size();
    const size_t cidx = lib.center();
    Emission<Abc> emission(lib.wlen(), weight_center, weight_decay, &sm);
    double oldloglike = loglike;
    loglike = 0.0;

    // Set progess bar workload to number of leapfrog steps times number of states
    if (prog_bar) prog_bar->Init(trainset.size());
    // Initialize sufficient statistics to pseudocount values
    Assign(priors, 0.0);
    for (size_t k = 0; k < probs.size(); ++k) {
      for (size_t j = 0; j < lib.wlen(); ++j) {
        for (size_t a = 0; a < Abc::kSize; ++a)
          probs[k][j][a] = pca * sm.p(a);
        probs[k][j][Abc::kAny] = 0.0;
      }
    }
    TransformToLog(lib);  // emission scores will be calculated in log-space

#pragma omp parallel for schedule(static)
    for (int n = 0; n < ntrain; ++n) {
      Vector<double> pp(lib.size(), 0.0);  // posterior P(z_n=k|c_n)
      double tmp = 0.0;

      // Compute posterior probs pp[k] of profile k for counts n
      tmp = CalculatePosteriorProbs(lib, emission, trainset[n], cidx, &pp[0]);
#pragma omp atomic
      loglike += tmp;

      // Update sufficient stasticics
      for (size_t k = 0; k < lib.size(); ++k) {
#pragma omp atomic
        priors[k] += pp[k];
        for (size_t j = 0; j < lib.wlen(); ++j) {
          for (size_t a = 0; a < Abc::kSize; ++a) {
#pragma omp atomic
            probs[k][j][a] += trainset[n].counts[j][a] * pp[k];
          }
        }
      }

      // Advance progress bar
      if (prog_bar) {
#pragma omp critical (advance_progress)
        prog_bar->Advance(1);
      }
    }
    loglike /= trainset.size() * emission.GetSumWeights();

    return loglike - oldloglike;
  }

  void MStep() {
    LOG(DEBUG) << StringifyRange(&priors[0], &priors[0] + priors.size());
    Normalize(&priors[0], priors.size());
    const size_t c = lib.center();
    for (size_t k = 0; k < lib.size(); ++k) {
      LOG(DEBUG) << probs[k];
      for (size_t j = 0; j < lib.wlen(); ++j)
        probs[k][j][Abc::kAny] = 1.0;
      Normalize(probs[k], 1.0);

      lib[k].is_log = false;
      lib[k].prior = priors[k];
      lib[k].probs = probs[k];
      std::copy(&probs[k][c][0], &probs[k][c][0] + Abc::kSizeAny, &lib[k].pc[0]);
    }
  }

  const TrainingSet& trainset;       // training set with counts we want to learn
  const SubstitutionMatrix<Abc>& sm; // substitution matrix needed for pseudocounts
  ContextLibrary<Abc> lib;           // current parameters set as context library
  Vector<double> priors;             // sufficient statistics for priors
  Vector<Profile<Abc> > probs;       // sufficient statistics for profile probs
  double weight_center;              // weight of central window position
  double weight_decay;               // exponential decay of window weights
  double pca;                        // pseudocount admix for profile probs
  double loglike;                    // current log-likelihood
};


// Driver function that runs EM-clustering until convergence, that is until
// absolute change in LL per column falls below threshold 'toll'. Return value
// is the number of EM iterations performed.
template<class Abc>
int RunClustering(EMClustering<Abc>& em,
                  double toll,
                  int max_iters,
                  FILE* fout = NULL) {
  scoped_ptr<ProgressBar> prog_bar;
  if (fout) {
    prog_bar.reset(new ProgressBar(fout, 30));
    fprintf(fout, "%-4s %-30s %9s %9s\n", "Iter", "E-Step", "Loglike", "+/-");
    fprintf(fout, "%s\n", std::string(55, '-').c_str());
  }

  int iter = 0;
  double delta = DBL_MAX;
  while ((iter == 0 || delta > toll) && iter < max_iters) {
    if (fout) { fprintf(fout, "%-4d ", iter); fflush(fout); }
    delta = em.EStep(prog_bar.get());
    if (fout) {
      if (iter == 0) fprintf(fout, " %9.6f\n", em.loglike);
      else fprintf(fout, " %9.6f %+9.6f\n", em.loglike, delta);
    }
    em.MStep();
    iter++;
  }

  return iter;
}

}  // namespace cs

#endif  // CS_EM_CLUSTERING_H_
