// Copyright 2009, Andreas Biegert

#ifndef CS_FUNC_H_
#define CS_FUNC_H_

#include "context_library-inl.h"
#include "crf-inl.h"
#include "emission.h"
#include "progress_bar.h"
#include "substitution_matrix-inl.h"
#include "training_sequence.h"
#include "training_profile.h"

namespace cs {

struct TrainingBlock {
    TrainingBlock() : beg(0), end(0), size(0), frac(0) {}

    TrainingBlock(size_t b, size_t e, size_t s, double f)
            : beg(b), end(e), size(s), frac(f) {}

    size_t beg;
    size_t end;
    size_t size;
    double frac;
};


template<class Abc, class TrainingPair>
struct ContextLibFunc {
    typedef std::vector<TrainingPair> TrainingSet;

    ContextLibFunc(const TrainingSet& tset,
                   const SubstitutionMatrix<Abc>& m,
                   double w_center = 1.6,
                   double w_decay = 0.85,
                   double t = 1.0)
            : trainset(tset),
              sm(m),
              weight_center(w_center),
              weight_decay(w_decay),
              tau(t) {}

    double operator() (const ContextLibrary<Abc>& lib,
                       ProgressBar* prog_bar = NULL) const {
        const int ntrain = trainset.size();
        const size_t cidx = lib.center();
        Emission<Abc> emission(lib.wlen(), weight_center, weight_decay, &sm);
        double loglike = 0.0;

        // Compute posterior probs p_zn[k] of profile k for counts n
#pragma omp parallel for schedule(static)
        for (int n = 0; n < ntrain; ++n) {
            Vector<double> pp(lib.size(), 0.0);  // posterior P(z_n=k|c_n)
            Vector<double> pa(Abc::kSize, 0.0);  // pseudocounts P(a|c_n)
            // First calcluate posteriors for each training window
            CalculatePosteriorProbs(lib, emission, trainset[n].x, cidx, &pp[0]);
            // Calculate pseudocounts p(a|c_n)
            for (size_t k = 0; k < lib.size(); ++k)
                for (size_t a = 0; a < Abc::kSize; ++a)
                    pa[a] += pp[k] * lib[k].pc[a];
            // Correct pseudocounts for admix factor 'tau'
            for (size_t a = 0; a < Abc::kSize; ++a)
                pa[a] = (1.0 - tau) * (trainset[n].x[cidx] == a ? 1.0 : 0.0) + tau * pa[a];
            // Calculate log-likelihood
            double loglike_n = 0.0;
            for (size_t a = 0; a < Abc::kSize; ++a)
                loglike_n += trainset[n].y[a] * (log(pa[a]) - log(sm.f(a)));
#pragma omp atomic
            loglike += loglike_n;
            // Advance progress bar
            if (prog_bar) {
#pragma omp critical (advance_progress)
                prog_bar->Advance(lib.size());
            }
        }

        return loglike;
    }

    const TrainingSet& trainset;
    const SubstitutionMatrix<Abc>& sm;
    double weight_center;
    double weight_decay;
    double tau;
};


template<class Abc, class TrainingPair>
struct CrfFunc {
    typedef std::vector<TrainingPair> TrainingSet;

    CrfFunc(const TrainingSet& tset, const SubstitutionMatrix<Abc>& m)
            : trainset(tset), sm(m) {}

    virtual ~CrfFunc() {}

    double operator() (const Crf<Abc>& crf, ProgressBar* prog_bar = NULL) const {
        double loglike = 0.0;
        const size_t center = crf.center();
        const int ntrain = trainset.size();

#pragma omp parallel for schedule(static)
        for (int n = 0; n < ntrain; ++n) {
            assert_eq(center, (trainset[n].x.length() - 1) / 2);
            Vector<double> pp(crf.size(), 0.0);  // posterior P(k|c_n)
            Vector<double> pa(Abc::kSize, 0.0);  // pseudocounts P(a|c_n)

            // Calculate posterior probability pp[k] of state k given count profile n
            double max = -DBL_MAX;
            for (size_t k = 0; k < crf.size(); ++k) {
                pp[k] = crf[k].bias_weight +
                    ContextScore(crf[k].context_weights, trainset[n].x, center, center);
                if (pp[k] > max) max = pp[k];  // needed for log-sum-exp trick
            }

            // Log-sum-exp trick begins here
            double sum = 0.0;
            for (size_t k = 0; k < crf.size(); ++k)
                sum += exp(pp[k] - max);

            double tmp = max + log(sum);
            for (size_t k = 0; k < crf.size(); ++k) {
                pp[k] = DBL_MIN + exp(pp[k] - tmp);
                // Calculate pseudocounts p(a|c_n)
                for (size_t a = 0; a < Abc::kSize; ++a)
                    pa[a] += crf[k].pc[a] * pp[k];
            }
            double loglike_n = 0.0;
            for (size_t a = 0; a < Abc::kSize; ++a) {
                pa[a] = MAX(DBL_MIN, pa[a]);
                loglike_n += trainset[n].y[a] * (log(pa[a]) - log(sm.p(a)));
            }
#pragma omp atomic
            loglike += loglike_n;

            // Advance progress bar
            if (prog_bar) {
#pragma omp critical (advance_progress)
                prog_bar->Advance(crf.size());
            }
        }

        return loglike;
    }

    const TrainingSet& trainset;
    const SubstitutionMatrix<Abc>& sm;
};


template<class Abc>
struct DerivCrfFuncIO {
    DerivCrfFuncIO(const Crf<Abc>& c)
            : crf(c),
              grad_loglike(crf.nweights(), 0.0),
              grad_prior(crf.nweights(), 0.0),
              loglike(-DBL_MAX),
              prior(-DBL_MAX) {}

    virtual ~DerivCrfFuncIO() {}

    Crf<Abc> crf;                 // CRF weights of current solution
    Vector<double> grad_loglike;  // partial derivatives of LL
    Vector<double> grad_prior;    // partial derivatives of prior
    double loglike;               // log of conditional likelihood
    double prior;                 // log of prior probability
};

template<class Abc>
struct DerivCrfFuncPrior {
    DerivCrfFuncPrior(double sb, double sc, double sd, double sp) : 
      sigma_bias(sb), sigma_context(sc), sigma_decay(sd), sigma_pc(sp) {}

    virtual double operator() (const Crf<Abc>& crf) const = 0; 

    virtual void CalculateGradient (
            const Crf<Abc>& crf,
            Vector<double>& grad,
            const TrainingBlock& block) const = 0;

    virtual ~DerivCrfFuncPrior() {}

    double sigma_bias;
    double sigma_context;
    double sigma_decay;
    double sigma_pc;
};

template<class Abc>
struct GaussianDerivCrfFuncPrior : public DerivCrfFuncPrior<Abc> {

    GaussianDerivCrfFuncPrior(
        double sb = 1.0, 
        double sc = 0.2, 
        double sd = 0.9,
        double sp = 10.0) : 
      DerivCrfFuncPrior<Abc>(sb, sc, sd, sp) {} 

    double operator() (const Crf<Abc>& crf) const {

        const double fac_bias = -0.5 / SQR(sigma_bias);
        // Precalculate factors for position specific sigmas in context window
        Vector<double> fac_context(crf.wlen());
        const int c = crf.center();
        for (size_t j = 0; j < crf.wlen(); ++j) {
            double tmp = sigma_context * pow(sigma_decay, fabs(static_cast<int>(j) - c));
            fac_context[j] = -0.5 / SQR(tmp);
        }
        const double fac_pc = -0.5 / SQR(sigma_pc);

        double prior = 0.0;
        for (size_t k = 0; k < crf.size(); ++k) {
            prior += fac_bias * SQR(crf[k].bias_weight);
            for(size_t j = 0; j < crf.wlen(); ++j) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    prior += fac_context[j] * SQR(crf[k].context_weights[j][a]);
            }
            const double* cw = crf[k].context_weights[c];
            for (size_t a = 0; a < Abc::kSize; ++a)
                prior += fac_pc * SQR(crf[k].pc_weights[a] - cw[a]);
        }
        return prior;
    }

    void CalculateGradient(
            const Crf<Abc>& crf,
            Vector<double>& grad,
            const TrainingBlock& block) const {

        const double fac_bias = -block.frac / SQR(sigma_bias);
        // Precalculate factors for position specific sigmas in context window
        Vector<double> fac_context(crf.wlen());
        const int c = crf.center();
        for (size_t j = 0; j < crf.wlen(); ++j) {
            double tmp = sigma_context * pow(sigma_decay, fabs(static_cast<int>(j) - c));
            fac_context[j] = -block.frac / SQR(tmp);
        }
        const double fac_pc = -block.frac / SQR(sigma_pc);

        Assign(grad, 0.0);  // reset gradient
        for (size_t k = 0, i = 0; k < crf.size(); ++k) {
            grad[i++] += fac_bias * crf[k].bias_weight;
            for(size_t j = 0; j < crf.wlen(); ++j) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    grad[i++] += fac_context[j] * crf[k].context_weights[j][a];
            }
            const double* cw = crf[k].context_weights[c];
            for (size_t a = 0; a < Abc::kSize; ++a)
                grad[i++] += fac_pc * (crf[k].pc_weights[a] - cw[a]);
        }
    }

    using DerivCrfFuncPrior<Abc>::sigma_bias;
    using DerivCrfFuncPrior<Abc>::sigma_context;
    using DerivCrfFuncPrior<Abc>::sigma_decay;
    using DerivCrfFuncPrior<Abc>::sigma_pc;
};

template<class Abc>
struct LassoDerivCrfFuncPrior : public DerivCrfFuncPrior<Abc> {

    LassoDerivCrfFuncPrior(
        double sb = 10.0, 
        double sc = 10.0, 
        double sd = 1.0,
        double sp = 10.0) : 
      DerivCrfFuncPrior<Abc>(sb, sc, sd, sp) {} 
   
    double operator() (const Crf<Abc>& crf) const {

        const double fac_bias = -1 / SQR(sigma_bias);
        // Precalculate factors for position specific sigmas in context window
        Vector<double> fac_context(crf.wlen());
        const int c = crf.center();
        for (size_t j = 0; j < crf.wlen(); ++j) {
            double tmp = sigma_context * pow(sigma_decay, fabs(static_cast<int>(j) - c));
            fac_context[j] = -1 / SQR(tmp);
        }
        const double fac_pc = -1 / SQR(sigma_pc);

        double prior = 0.0;
        for (size_t k = 0; k < crf.size(); ++k) {
            prior += fac_bias * fabs(crf[k].bias_weight);
            for(size_t j = 0; j < crf.wlen(); ++j) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    prior += fac_context[j] * fabs(crf[k].context_weights[j][a]);
            }
            const double* cw = crf[k].context_weights[c];
            for (size_t a = 0; a < Abc::kSize; ++a)
                prior += fac_pc * fabs(crf[k].pc_weights[a] - cw[a]);
        }
        return prior;
    }

    void CalculateGradient(
            const Crf<Abc>& crf,
            Vector<double>& grad,
            const TrainingBlock& block) const {

        const double fac_bias = -block.frac / SQR(sigma_bias);
        // Precalculate factors for position specific sigmas in context window
        Vector<double> fac_context(crf.wlen());
        const int c = crf.center();
        for (size_t j = 0; j < crf.wlen(); ++j) {
            double tmp = sigma_context * pow(sigma_decay, fabs(static_cast<int>(j) - c));
            fac_context[j] = -block.frac / SQR(tmp);
        }
        const double fac_pc = -block.frac / SQR(sigma_pc);

        Assign(grad, 0.0);  // reset gradient
        for (size_t k = 0, i = 0; k < crf.size(); ++k) {
            grad[i++] += fac_bias * SIGN(crf[k].bias_weight);
            for(size_t j = 0; j < crf.wlen(); ++j) {
                for (size_t a = 0; a < Abc::kSize; ++a)
                    grad[i++] += fac_context[j] * SIGN(crf[k].context_weights[j][a]);
            }
            const double* cw = crf[k].context_weights[c];
            for (size_t a = 0; a < Abc::kSize; ++a)
                grad[i++] += fac_pc * SIGN(crf[k].pc_weights[a] - cw[a]);
        }
    }

    using DerivCrfFuncPrior<Abc>::sigma_bias;
    using DerivCrfFuncPrior<Abc>::sigma_context;
    using DerivCrfFuncPrior<Abc>::sigma_decay;
    using DerivCrfFuncPrior<Abc>::sigma_pc;
};

template<class Abc, class TrainingPair>
struct DerivCrfFunc : public CrfFunc<Abc, TrainingPair> {
    typedef std::vector<TrainingPair> TrainingSet;

    DerivCrfFunc(const TrainingSet& trainset,
                 const SubstitutionMatrix<Abc>& m,
                 DerivCrfFuncPrior<Abc>& p)
            : CrfFunc<Abc, TrainingPair>(trainset, m),
              shuffle(trainset.size()),
              prior(p) {
        for (size_t i = 0; i < trainset.size(); ++i) shuffle[i] = i;
    }

    TrainingBlock GetBlock(size_t b, size_t nblocks) const {
        assert(b < nblocks);
        size_t block_size = trainset.size() / nblocks;
        size_t beg = b * block_size;
        size_t end = (b == nblocks - 1) ? trainset.size() : (b + 1) * block_size;
        size_t size = end - beg;
        double frac = static_cast<double>(size) / trainset.size();
        return TrainingBlock(beg, end, size, frac);
    }

    void df(DerivCrfFuncIO<Abc>& s,
            size_t b = 0,  // index of training block
            size_t nblocks = 1,  // total number of training blocks
            ProgressBar* prog_bar = NULL) const {
        assert(b < nblocks);
        const TrainingBlock block(GetBlock(b, nblocks));
        const int n_beg = block.beg;
        const int n_end = block.end;
        const size_t center = s.crf.center();
        Matrix<double> mpp(block.size, s.crf.size(), 0.0);  // posterior P(k|c_n)
        Matrix<double> mpa(block.size, Abc::kSize, 0.0);    // pseudocounts P(a|c_n)
        double loglike = 0.0;

#pragma omp parallel for schedule(static)
        for (int n = n_beg; n < n_end; ++n) {
            const TrainingPair& tpair = trainset[shuffle[n]];
            //assert_eq(center, (tpair.x.length() - 1) / 2);
            double* pp = &mpp[n - block.beg][0];
            double* pa = &mpa[n - block.beg][0];

            // Calculate posterior probability pp[k] of state k given count profile n
            double max = -DBL_MAX;
            for (size_t k = 0; k < s.crf.size(); ++k) {
                pp[k] = s.crf[k].bias_weight +
                    ContextScore(s.crf[k].context_weights, tpair.x, center, center);
                if (pp[k] > max) max = pp[k];  // needed for log-sum-exp trick
            }

            // Log-sum-exp trick begins here
            double sum = 0.0;
            for (size_t k = 0; k < s.crf.size(); ++k)
                sum += exp(pp[k] - max);
            double tmp = max + log(sum);

            for (size_t k = 0; k < s.crf.size(); ++k) {
                pp[k] = DBL_MIN + exp(pp[k] - tmp);
                // Calculate pseudocounts p(a|c_n)
                for (size_t a = 0; a < Abc::kSize; ++a)
                    pa[a] += s.crf[k].pc[a] * pp[k];
            }

            long double loglike_n = 0.0;
            for (size_t a = 0; a < Abc::kSize; ++a) {
                pa[a] = MAX(DBL_MIN, pa[a]);
                loglike_n += tpair.y[a] * (log(pa[a]) - log(sm.p(a)));
            }
#pragma omp atomic
            loglike += loglike_n;
        }

        s.loglike += loglike;
        s.prior   += block.frac * prior(s.crf);

        CalculateLikelihoodGradient(trainset, block, s.crf, mpp, mpa, s.grad_loglike, prog_bar);
        prior.CalculateGradient(s.crf, s.grad_prior, block);
    }

    // This is the performance critical method in HMC sampling. It accounts for about
    // 90% of the runtime in profiling.
    void CalculateLikelihoodGradient(const std::vector<TrainingSequence<Abc> >& tset, 
                                     const TrainingBlock& block,
                                     const Crf<Abc>& crf,
                                     const Matrix<double>& mpp,
                                     const Matrix<double>& mpa,
                                     Vector<double>& grad,
                                     ProgressBar* prog_bar) const {
        const size_t wlen = crf.wlen();
        const int nstates = crf.size();
        Assign(grad,  0.0);  // reset gradient

        // We perform parallelization over training points instead of CRF states because
        // this way we don't have to protect access to the gradient vector.
#pragma omp parallel for schedule(static)
        for (int k = 0; k < nstates; ++k) {
            const double* pc = &(crf[k].pc[0]);
            double fit, sum = 0.0;
            size_t offset = k * (1 + (wlen + 1) * Abc::kSize);

            for (size_t n = block.beg; n < block.end; ++n) {
                const size_t m = n - block.beg;  // index of training point 'n' in matrices
                const double* pa = &mpa[m][0];   // predicted pseudocounts
                const TrainingSequence<Abc>& tseq = tset[shuffle[n]];
                size_t i = offset;

                // Precumpute fit
                fit = 0.0;
                for (size_t a = 0; a < Abc::kSize; ++a)
                    fit += tseq.y[a] * (pc[a] / pa[a] - 1.0);
                double mpp_fit = mpp[m][k] * fit;
                mpp_fit = MIN(DBL_MAX, mpp_fit);

                // Update gradient of bias weight
                grad[i++] += mpp_fit;

                // Update gradient terms of context weights
                for(size_t j = 0; j < wlen; ++j) {
                    if (tseq.x[j] != Abc::kAny)
                        grad[i + tseq.x[j]] += mpp_fit;
                    i += Abc::kSize;
                }

                // Precompute sum needed for gradient update of pseudocounts weights
                sum = 0.0;
                for (size_t a = 0; a < Abc::kSize; ++a)
                    sum += pc[a] * tseq.y[a] / pa[a];
                sum = MIN(DBL_MAX, sum);


                // Update gradient terms of pseudocount weights
                for (size_t a = 0; a < Abc::kSize; ++a) {
                    double d = tseq.y[a] / pa[a];
                    d = MIN(DBL_MAX, d);
                    grad[i++] += mpp[m][k] * pc[a] * (d - sum);
                }
            }
            // Advance progress bar
            if (prog_bar) {
#pragma omp critical (advance_progress)
                prog_bar->Advance(block.end - block.beg);
            }
        }
    }

    // This is the performance critical method in HMC sampling. It accounts for about
    // 90% of the runtime in profiling.
    void CalculateLikelihoodGradient(const std::vector<TrainingProfile<Abc> >& tset, 
                                     const TrainingBlock& block,
                                     const Crf<Abc>& crf,
                                     const Matrix<double>& mpp,
                                     const Matrix<double>& mpa,
                                     Vector<double>& grad,
                                     ProgressBar* prog_bar) const {
        const size_t wlen = crf.wlen();
        const int nstates = crf.size();
        Assign(grad,  0.0);  // reset gradient

        // We perform parallelization over training points instead of CRF states because
        // this way we don't have to protect access to the gradient vector.
#pragma omp parallel for schedule(static)
        for (int k = 0; k < nstates; ++k) {
            const double* pc = &(crf[k].pc[0]);
            double fit, sum = 0.0;
            size_t offset = k * (1 + (wlen + 1) * Abc::kSize);

            for (size_t n = block.beg; n < block.end; ++n) {
                const size_t m = n - block.beg;  // index of training point 'n' in matrices
                const double* pa = &mpa[m][0];   // predicted pseudocounts
                const TrainingProfile<Abc>& tprof = tset[shuffle[n]];
                size_t i = offset;

                // Precumpute fit
                fit = 0.0;
                for (size_t a = 0; a < Abc::kSize; ++a)
                    fit += tprof.y[a] * (pc[a] / pa[a] - 1.0);
                double mpp_fit = mpp[m][k] * fit;
                mpp_fit = MIN(DBL_MAX, mpp_fit);

                // Update gradient of bias weight
                grad[i++] += mpp_fit;

                // Update gradient terms of context weights
                for(size_t j = 0; j < wlen; ++j)
                    for (size_t a = 0; a < Abc::kSize; ++a)
                        grad[i++] += tprof.x.counts[j][a] * mpp_fit;						

                // Precompute sum needed for gradient update of pseudocounts weights
                sum = 0.0;
                for (size_t a = 0; a < Abc::kSize; ++a)
                    sum += pc[a] * tprof.y[a] / pa[a];
                sum = MIN(DBL_MAX, sum);

                // Update gradient terms of pseudocount weights
                for (size_t a = 0; a < Abc::kSize; ++a) {
                    double d = tprof.y[a] / pa[a];
                    d = MIN(DBL_MAX, d);
                    grad[i++] += mpp[m][k] * pc[a] * (d - sum);
                }
            }
            // Advance progress bar
            if (prog_bar) {
#pragma omp critical (advance_progress)
                prog_bar->Advance(block.end - block.beg);
            }
        }


    }

    using CrfFunc<Abc, TrainingPair>::trainset;
    using CrfFunc<Abc, TrainingPair>::sm;
    std::vector<int> shuffle;
    DerivCrfFuncPrior<Abc>& prior;
}; 



}  // namespace cs

#endif  // CS_FUNC_H_
