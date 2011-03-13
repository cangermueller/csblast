// Copyright 2009, Andreas Biegert

#ifndef CS_SGD_H_
#define CS_SGD_H_

#include "crf-inl.h"
#include "func.h"
#include "progress_bar.h"

namespace cs {

// Parameter wrapper for various SGD parameters
struct SgdParams {
    SgdParams()
            : nblocks(100),
              eta(0.001),
              mu(0.01),
              rho(0.5),
              max_eta(1.0),
              gamma(0.9),
              toll(1e-3),
              min_epochs(10),
              max_epochs(100),
              seed(0) {}

    size_t nblocks;     // number of training blocks
    double eta;         // learning rate eta for context and pseudocount weights
    double mu;          // meta learning rate
    double rho;         // lower bound for multiplier in learning rate adaption
    double max_eta;     // upper bound for learning rates
    double gamma;       // parameter governing exponential average of derivatives
    double toll;        // LL change for convergence
    int min_epochs;     // minimal number of epochs
    int max_epochs;     // maximal number of epochs
    unsigned int seed;  // seed for rng that shuffles training set after each epoch
};


template<class Abc>
struct SgdState : public DerivCrfFuncIO<Abc> {
    SgdState(const Crf<Abc>& c)
            : DerivCrfFuncIO<Abc>(c),
              steps(0),
              eta(c.nweights(), 0.0),
              avg(c.nweights(), 0.0),
              grad_prev(c.nweights(), 0.0) {}

    size_t steps;              // number of SGD steps already performed
    Vector<double> eta;        // learning rates eta for each CRF weight
    Vector<double> avg;        // exponential average of learning rates
    Vector<double> grad_prev;  // previous gradient of likelihood and prior combined
};


template<class Abc, class TrainingPair>
struct Sgd {
    Sgd(const DerivCrfFunc<Abc, TrainingPair>& tf,
        const SgdParams& params)
            : func(tf),
              nblocks(params.nblocks),
              eta(params.eta),
              mu(params.mu),
              rho(params.rho),
              max_eta(params.max_eta),
              gamma(params.gamma),
              ran(params.seed) {}

    // Shuffles training set and then runs one epoche of stochastic gradient descent
    // comprising of 'nblocks' individual gradient descent steps.
    void operator() (SgdState<Abc>& s, ProgressBar* prog_bar = NULL) {
        s.loglike = 0.0;
        s.prior = 0.0;
        // Shuffle training set before each epoch
        random_shuffle(func.shuffle.begin(), func.shuffle.end(), ran);

        for (size_t b = 0; b < nblocks; ++b) {
            // Save previous gradient of likelihood and prior as combined vector
            for (size_t i = 0; i < s.grad_prev.size(); ++i)
                s.grad_prev[i] = s.grad_loglike[i] + s.grad_prior[i];
            // Calculate gradient and increment likelihood based on training block 'b'
            func.df(s, b, nblocks, prog_bar);
            // Udpate averages of partial derivatives
            UpdateAverages(s);
            // Udpate learning rates
            UpdateLearningRates(s);
            // Updates CRF weights based on learning rates and gradient
            UpdateCRF(s);

            s.steps++;
        }
    }

    // Moves CRF weights along the gradient direction with parameter-specific
    // step sizes given in 'eta' vector.
    void UpdateCRF(SgdState<Abc>& s) {
        for (size_t k = 0, i = 0; k < s.crf.size(); ++k) {
            s.crf[k].bias_weight += s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
            ++i;
            for (size_t j = 0; j < s.crf.wlen(); ++j)
                for (size_t a = 0; a < Abc::kSize; ++a, ++i)
                    s.crf[k].context_weights[j][a] += s.eta[i] * (s.grad_loglike[i] +
                                                                  s.grad_prior[i]);
            for (size_t a = 0; a < Abc::kSize; ++a, ++i)
                s.crf[k].pc_weights[a] += s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
            UpdatePseudocounts(s.crf[k]);
        }
    }

    // Updates all learning rates based on the current gradient, the previous
    // gradient and the exponential average vector.
    void UpdateLearningRates(SgdState<Abc>& s) {
        if (s.steps == 0) {
            Assign(s.eta, eta); // start with constant eta
        } else {
            double tmp = 0.0;
            for (size_t i = 0; i < s.eta.size(); ++i) {
                tmp = (s.grad_loglike[i] + s.grad_prior[i]) * s.grad_prev[i] / s.avg[i];
                s.eta[i] = MIN(max_eta, s.eta[i] * MAX(rho, 1.0 + mu * tmp));
            }
        }
    }

    // Updates exponential averages with new gradient
    void UpdateAverages(SgdState<Abc>& s) {
        double tmp = s.steps == 0 ? 1.0 : 1.0 - gamma;
        for (size_t i = 0; i < s.avg.size(); ++i)
            s.avg[i] = gamma * s.avg[i] + tmp * SQR(s.grad_loglike[i] + s.grad_prior[i]);
    }

    DerivCrfFunc<Abc, TrainingPair> func;  // training set function
    size_t nblocks;     // number of training blocks
    double eta;         // learning rate eta for context and pseudocount weights
    double mu;          // meta learning rate
    double rho;         // lower bound for multiplier in learning rate adaption
    double max_eta;     // upper bound for learning rates
    double gamma;       // parameter governing exponential average of derivatives
    Ran ran;            // RNG for shuffling of training set
};


template<class Abc, class TrainingPair>
struct SgdOptimizer {
    SgdOptimizer(const DerivCrfFunc<Abc, TrainingPair>& tf,
                 const CrfFunc<Abc, TrainingPair>& vf,
                 const SgdParams& params)
            : sgd(tf, params),
              func(vf),
              toll(params.toll),
              min_epochs(params.min_epochs),
              max_epochs(params.max_epochs) {}

    double Optimize(Crf<Abc>& crf, FILE* fout = NULL) {
        scoped_ptr<ProgressBar> prog_bar;
        SgdState<Abc> s(crf);
        int epoch = 0, nconv = 0;
        double best_loglike = -DBL_MAX, delta = DBL_MAX, old_loglike, val_loglike;

        if (fout) {
            prog_bar.reset(new ProgressBar(fout, 16));
            fprintf(fout, "%-5s %-16s  %8s %7s  %9s  %7s\n", "Epoch",
                    "Gradient descent", "LL-Train", "+/-", "Prior", "LL-Val");
            fprintf(fout, "%s\n", std::string(60, '-').c_str());
        }

        while ((nconv < kMaxConvBumps || epoch < min_epochs) && epoch < max_epochs) {
            // Print first part of table row
            if (fout) {
                fprintf(fout, "%-4d  ", epoch + 1); fflush(fout);
                prog_bar->Init(sgd.func.trainset.size() * crf.size());
            }
            // Save last likelihood for calculation of delta
            old_loglike = s.loglike;
            // Run on epoche of SGD
            sgd(s, prog_bar.get());
            // Normalize likelihood and prior to user friendly scale
            s.loglike /= sgd.func.trainset.size();
            s.prior /= crf.nweights();
            delta = s.loglike - old_loglike;
            // Keep track of how many times we were under convergence threshold
            if (fabs(delta) > toll) nconv = 0;
            else ++nconv;
            // Calculate likelihood on validation set
            val_loglike = func(s.crf) / func.trainset.size();
            if (val_loglike > best_loglike) {
                best_loglike = val_loglike;
                crf = s.crf;
            }
            // Print second part of table row
            if (fout) {
                fprintf(fout, " %8.4f %+8.4f %10.4f %8.4f\n",
                        s.loglike, epoch == 0 ? 0.0 : delta, s.prior, val_loglike);
            }
            epoch++;
        }

        return best_loglike;
    }

    static const int kMaxConvBumps = 3;

    Sgd<Abc, TrainingPair> sgd;       // SGD algorithm encapsulation
    CrfFunc<Abc, TrainingPair> func;  // validation set function
    double toll;                      // LL change for convergence
    int min_epochs;                   // minimal number of epochs
    int max_epochs;                   // maximal number of epochs
};


}  // namespace cs

#endif  // CS_SGD_H_
