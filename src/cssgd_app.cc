// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf-inl.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "sgd.h"
#include "matrix_pseudocounts-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"
#include "training_profile.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

// CSSgdAppOptions
struct CSSgdAppOptions {
    CSSgdAppOptions() { Init(); }

    void Init() {
        nstates        = 100;
        blosum_type    = "BLOSUM62";
        pc_init        = 0.75f;
        gauss_init     = 0.0;
        prior          = 2;
        weight_center  = 1.6;
        weight_decay   = 0.85;
        neff_ext       = "seq";
        neff_nsamples  = 100;
        neff_pc        = 1.0;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (trainfile.empty()) throw Exception("No training set provided!");
        if (valfile.empty()) throw Exception("No validation set provided!");
        if (crffile_vset.empty()) throw Exception("No output file for CRF provided!");
        if (prior < 1 || prior > 3) throw Exception("Invalid prior!");
        if (pc_init <= 0 || pc_init > 1.0) throw Exception("Pseudocounts admix for initialization invalid!");
        if (neff_pc <= 0 || neff_pc > 1.0) throw Exception("Pseudocounts admix for computing the Neff invalid!");
        if (sgd.eta_mode < 1 || sgd.eta_mode > 2) throw Exception("Invalid mode for updating the learning rate eta!");
        if (sgd.eta_decay < 1) throw Exception("Eta decay must greater/equal one!");
    }

    void PrintOptions(FILE* out) const {
        fprintf(out, "  %3s %-25s: %s\n", "-t,", "--trainset", trainfile.c_str()); 
        fprintf(out, "  %3s %-25s: %s\n", "-v,", "--valset", valfile.c_str()); 
        fprintf(out, "  %3s %-25s: %s\n", "-o,", "--outfile", crffile_vset.c_str()); 
        fprintf(out, "  %3s %-25s: %s\n", "-O,", "--outfile-tset", crffile_tset.c_str()); 
        fprintf(out, "  %3s %-25s: %s\n", "-R,", "--progress", outfile.c_str()); 
        fprintf(out, "  %3s %-25s: %zu\n", "-K,", "--states", nstates); 
        fprintf(out, "  %3s %-25s: %zu\n", "-N,", "--epochs", sgd.max_epochs); 
        fprintf(out, "  %3s %-25s: %.3g\n", "-t,", "--toll", sgd.toll); 
        fprintf(out, "  %3s %-25s: %.3g\n", "-T,", "--early-delta", sgd.early_delta); 
        fprintf(out, "  %3s %-25s: %.3g\n", "", "--min-ll", sgd.min_ll); 
        fprintf(out, "  %3s %-25s: %zu\n", "", "--min-ll-repeats", sgd.min_ll_repeats); 
        fprintf(out, "\n");

        fprintf(out, "  %3s %-25s: %zu\n", "-E,", "--eta-mode", sgd.eta_mode); 
        fprintf(out, "  %3s %-25s: %.3g\n", "-e,", "--eta", sgd.eta_init); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-D,", "--eta-decay", sgd.eta_decay); 
        fprintf(out, "  %3s %-25s: %.3g\n", "", "--eta-reinit", sgd.eta_reinit); 
        fprintf(out, "  %3s %-25s: %zu\n", "", "--eta-reinit-num", sgd.eta_reinit_num); 
        fprintf(out, "  %3s %-25s: %.3g\n", "", "--eta-reinit-delta", sgd.eta_reinit_delta); 
        fprintf(out, "  %3s %-25s: %zu\n", "-B,", "--blocks", sgd.nblocks); 
        fprintf(out, "  %3s %-25s: %zu\n", "-P,", "--prior", prior); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-b,", "--sigma-bias", sgd.sigma_bias); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-c,", "--sigma-context", sgd.sigma_context); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-C,", "--sigma-context-pos", sgd.sigma_context_pos_max); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-d,", "--sigma-decay", sgd.sigma_decay); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-p,", "--sigma-pc", sgd.sigma_pc_max); 
        fprintf(out, "  %3s %-25s: %zu\n", "-q,", "--sigma-relax-epoch", sgd.sigma_relax_epoch); 
        fprintf(out, "  %3s %-25s: %zu\n", "", "--sigma-relax-steps", sgd.sigma_relax_steps); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-s,", "--context-penalty", sgd.context_penalty); 
        fprintf(out, "  %3s %-25s: %zu\n", "-S,", "--context-penalty-epoch", sgd.context_penalty_epoch); 
        fprintf(out, "  %3s %-25s: %zu\n", "" , "--context-penalty-steps", sgd.context_penalty_steps); 
        fprintf(out, "\n");

        fprintf(out, "  %3s %-25s: %s\n", "-m,", "--model", modelfile.c_str()); 
        fprintf(out, "  %3s %-25s: %.2f\n", "", "--weight-center", weight_center); 
        fprintf(out, "  %3s %-25s: %.2f\n", "", "--weight-decay", weight_decay); 
        fprintf(out, "  %3s %-25s: %.2f\n", "-g,", "--gauss-init", gauss_init); 
        fprintf(out, "  %3s %-25s: %.2f\n", "", "--pc-init", pc_init); 
        fprintf(out, "  %3s %-25s: %u\n", "-r,", "--seed", sgd.seed); 
        fprintf(out, "\n");

        fprintf(out, "  %3s %-25s: %s\n", "", "--neff-dir", neff_dir.c_str()); 
        fprintf(out, "  %3s %-25s: %s\n", "", "--neff-ext", neff_ext.c_str()); 
        fprintf(out, "  %3s %-25s: %zu\n", "", "--neff-nsamples", neff_nsamples); 
        fprintf(out, "  %3s %-25s: %.2f\n", "", "--neff-pc", neff_pc); 
    }


    // Input file with training set.
    string trainfile;
    // Input file with validation.
    string valfile;
    // Output file for the best CRF on the validation set
    string crffile_vset;
    // Output file for the last CRF on the training set
    string crffile_tset;
    // Output file for progress table.
    string outfile;
    // Crf input file with context library or CRF for initialization.
    string modelfile;
    // The number of states in the HMM to train.
    size_t nstates;
    // BLOSUM matrix for pseudocount generation.
    string blosum_type;
    // Constant admix for CRF initialization
    double pc_init;
    // Sigma for CRF init from gaussian
    double gauss_init;
    // Prior to be used for calculating the likelihood
    size_t prior;
    // Weight of central column in multinomial emission
    double weight_center;
    // Exponential decay of window weights
    double weight_decay;
    // Wrapper for SGD parameters
    SgdParams sgd;
    // Directory with sequences or profiles to be used for calculating the Neff
    string neff_dir;
    // File extension of sequences or profiles to be used for calculating the Neff
    string neff_ext;
    // Number of samples for calculating the Neff
    size_t neff_nsamples;
    // PC admixture for calculating the Neff
    double neff_pc;

};  // CSSgdAppOptions









// CSSgdRunner
template<class Abc, class TrainingPairT, class TrainingPairV>
struct CSSgdRunner {
    typedef vector<TrainingPairT> TrainingSetT;
    typedef vector<TrainingPairV> TrainingSetV;

    CSSgdRunner(FILE* const out, const CSSgdAppOptions& opts) : 
        out_(out),
        opts_(opts) {}


    // Initializes substitution matrix (specialized by alphabet type).
    void InitSubstitutionMatrix() {
        BlosumType type = BlosumTypeFromString(opts_.blosum_type);
        sm_.reset(new BlosumMatrix(type));
    }

    // Reads training and validation set.
    void ReadTrainingData() {
        FILE* fin;
        fprintf(out_, "Reading training set from %s ...\n",
                GetBasename(opts_.trainfile).c_str());
        fin = fopen(opts_.trainfile.c_str(), "r");
        if (!fin)
            throw Exception("Can't read training set from '%s'!", opts_.trainfile.c_str());
        ReadAll(fin, trainset_);
        fclose(fin);
        fprintf(out_, "%zu records read.\n", trainset_.size());
        if (trainset_.size() == 0)
            throw Exception("Training set empty!");

        fprintf(out_, "Reading validation set from %s ...\n",
                GetBasename(opts_.valfile).c_str());
        fin = fopen(opts_.valfile.c_str(), "r");
        if (!fin)
            throw Exception("Can't read validation set from '%s'!",
                            opts_.valfile.c_str());
        ReadAll(fin, valset_);
        fclose(fin);
        fprintf(out_, "%zu records read.\n", valset_.size());
        if (valset_.size() == 0)
            throw Exception("Validation set empty!");
        size_t s = MIN(trainset_.size(), valset_.size());
        if (s < opts_.sgd.nblocks) {
          opts_.sgd.nblocks = s;
          fprintf(out_, "Warning: block size changed to %zu blocks!\n", opts_.sgd.nblocks);
        }
    }

    // Initializes CRF from context library, training data, or jumpstart file.
    void InitCrf() {
        // Read CRF from jumpstart file?
        if (!opts_.modelfile.empty() && GetFileExt(opts_.modelfile) == "crf") {
            fprintf(out_, "Reading CRF from %s ...\n", GetBasename(opts_.modelfile).c_str());

            FILE* fin = fopen(opts_.modelfile.c_str(), "r");
            if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
            crf_.reset(new Crf<Abc>(fin));
            fclose(fin);

        } else if (!opts_.modelfile.empty()) {  // init from profile library
            fprintf(out_, "Initializing CRF with context library read from %s ...\n",
                    GetBasename(opts_.modelfile).c_str());

            FILE* fin = fopen(opts_.modelfile.c_str(), "r");
            if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
            ContextLibrary<Abc> lib(fin);
            fclose(fin);
            size_t nstates = opts_.nstates == 0 ? lib.size() : opts_.nstates;
            LibraryBasedCrfInit<Abc> init(lib, opts_.weight_center, opts_.weight_decay);
            crf_.reset(new Crf<Abc>(nstates, lib.wlen(), init));

        } else if (opts_.gauss_init != 0.0) {  // init by sampling weights from gaussian
            fputs("Initializing CRF by sampling weights from gaussian ...\n", out_);
            GaussianCrfInit<Abc> init(opts_.gauss_init, *sm_, opts_.sgd.seed);
            crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));

        } else {  // init by sampling from training set
            fputs("Initializing CRF by sampling from training set ...\n", out_);
            MatrixPseudocounts<Abc> pc(*sm_);
            ConstantAdmix admix(opts_.pc_init);
            SamplingCrfInit<Abc, TrainingPairT> init(trainset_, pc, admix, *sm_, 
                opts_.sgd.seed, opts_.weight_center, opts_.weight_decay);
            crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));
        }
    }

    // Initializes the samples for computing the Neff
    void InitNeffSamples() {
        fprintf(out_, "Initializing samples for calculating the Neff...\n");
        Ran ran(opts_.sgd.seed);
        if (!opts_.neff_dir.empty()) {
            vector<string> files;
            GetAllFiles(opts_.neff_dir, files, opts_.neff_ext);
            random_shuffle(files.begin(), files.end(), ran);
            size_t nsamples = MIN(opts_.neff_nsamples, files.size());
            bool profiles = opts_.neff_ext == "prf";
            for (size_t i = 0; i < nsamples; ++i) {
                const char* fn = PathCat(opts_.neff_dir, files[i]);                
                FILE* fin = fopen(fn, "r");
                if (!fin) throw Exception("Unable to read file '%s'!", fn);
                if (profiles)
                    neff_samples_.push_back(CountProfile<Abc>(fin));
                else 
                    neff_samples_.push_back(CountProfile<Abc>(Sequence<Abc>(fin)));
            }
        }
        if (neff_samples_.size() == 0) {
            size_t wlen = trainset_.front().x.length();
            for (size_t i = 0; i < opts_.neff_nsamples; ++i) {
                CountProfile<Abc> sample(wlen * kNeffWindowsPerSample);
                for (size_t j = 0; j < sample.length(); j += wlen) {
                    const CountProfile<Abc> cp(valset_[ran(valset_.size())].x);
                    sample.counts.Insert(j, cp.counts);
                }
                neff_samples_.push_back(sample);
            }             
        }
        fprintf(out_, "%zu samples initialized.\n", neff_samples_.size());
    }
     

    // Runs stochastic gradient descent.
    int Run() {
        fputs("Running stochastic gradient descent with following parameters:\n", out_);
        opts_.PrintOptions(out_);
        fputs("\n", out_);

        InitSubstitutionMatrix();
        ReadTrainingData();
        InitCrf();
        InitNeffSamples();

        FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");
        fputs("\n", fout);

        scoped_ptr<DerivCrfFuncPrior<Abc> > prior;
        if (opts_.prior == 1)
            prior.reset(new LassoDerivCrfFuncPrior<Abc>(
                opts_.sgd.sigma_bias,
                opts_.sgd.sigma_context, 
                opts_.sgd.sigma_decay, 
                opts_.sgd.sigma_pc_max, 
                opts_.sgd.context_penalty));
        else if (opts_.prior == 2) {
            prior.reset(new GaussianDerivCrfFuncPrior<Abc>(
                opts_.sgd.sigma_bias,
                opts_.sgd.sigma_context, 
                opts_.sgd.sigma_decay, 
                opts_.sgd.sigma_pc_max,
                opts_.sgd.context_penalty));
        } else {
            prior.reset(new UnsymmetricDerivCrfFuncPrior<Abc>(
                opts_.sgd.sigma_bias,
                opts_.sgd.sigma_context, 
                opts_.sgd.sigma_context_pos_max, 
                opts_.sgd.sigma_decay, 
                opts_.sgd.sigma_pc_max,
                opts_.sgd.context_penalty));
        }

        CrfFunc<Abc, TrainingPairV> val_func(valset_, *sm_);
        DerivCrfFunc<Abc, TrainingPairT> train_func(trainset_, *sm_, *prior);
        SgdOptimizer<Abc, TrainingPairT, TrainingPairV> sgd(train_func, val_func, 
            opts_.sgd, neff_samples_, opts_.neff_pc);
        sgd.crffile_tset = opts_.crffile_tset;
        sgd.crffile_vset = opts_.crffile_vset;
        sgd.Optimize(*crf_, fout);

        if (!opts_.outfile.empty()) fclose(fout);

        fprintf(out_, "\nWrote best CRF on validation set to %s!\n", opts_.crffile_vset.c_str());
        if (!opts_.crffile_tset.empty())
          fprintf(out_, "Wrote last CRF on training set to %s!\n", opts_.crffile_tset.c_str());

        return 0;
    }

    FILE* const out_;
    CSSgdAppOptions opts_;
    TrainingSetT trainset_;
    TrainingSetV valset_;
    scoped_ptr<Crf<Abc> > crf_;
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;
    vector<CountProfile<Abc> > neff_samples_;

    static const size_t kNeffWindowsPerSample = 3;

};  //CSSgdRunner









// CSSgdApp
template<class Abc>
class CSSgdApp : public Application {
  protected:
    // Options
    CSSgdAppOptions opts_;

    // Runs the csbuild application.
    virtual int Run();
    // Parses command line options.
    virtual void ParseOptions(GetOpt_pp& ops);
    // Prints options summary to stream.
    virtual void PrintOptions() const;
    // Prints short application description.
    virtual void PrintBanner() const;
    // Prints usage banner to stream.
    virtual void PrintUsage() const;

};  // CSSgdApp

template<class Abc>
void CSSgdApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "trainset", opts_.trainfile, opts_.trainfile);
    ops >> Option('j', "valset", opts_.valfile, opts_.valfile);
    ops >> Option('o', "outfile", opts_.crffile_vset, opts_.crffile_vset);
    ops >> Option('O', "outfile-tset", opts_.crffile_tset, opts_.crffile_tset);
    ops >> Option('R', "progress", opts_.outfile, opts_.outfile);
    ops >> Option('K', "states", opts_.nstates, opts_.nstates);
    ops >> Option('N', "epochs", opts_.sgd.max_epochs, opts_.sgd.max_epochs);
    ops >> Option('t', "toll", opts_.sgd.toll, opts_.sgd.toll);
    ops >> Option('T', "early-delta", opts_.sgd.early_delta, opts_.sgd.early_delta);
    ops >> Option(' ', "min-ll", opts_.sgd.min_ll, opts_.sgd.min_ll);
    ops >> Option(' ', "min-ll-repeats", opts_.sgd.min_ll_repeats, opts_.sgd.min_ll_repeats);

    ops >> Option('E', "eta-mode", opts_.sgd.eta_mode, opts_.sgd.eta_mode);
    ops >> Option('e', "eta", opts_.sgd.eta_init, opts_.sgd.eta_init);
    ops >> Option('D', "eta-decay", opts_.sgd.eta_decay, opts_.sgd.eta_decay);
    ops >> Option(' ', "eta-reinit", opts_.sgd.eta_reinit, opts_.sgd.eta_reinit);
    ops >> Option(' ', "eta-reinit-num", opts_.sgd.eta_reinit_num, opts_.sgd.eta_reinit_num);
    ops >> Option(' ', "eta-reinit-delta", opts_.sgd.eta_reinit_delta, opts_.sgd.eta_reinit_delta);
    ops >> Option('B', "blocks", opts_.sgd.nblocks, opts_.sgd.nblocks);
    ops >> Option('P', "prior", opts_.prior, opts_.prior);
    ops >> Option('b', "sigma-bias", opts_.sgd.sigma_bias, opts_.sgd.sigma_bias);
    ops >> Option('c', "sigma-context", opts_.sgd.sigma_context, opts_.sgd.sigma_context);
    ops >> Option('C', "sigma-context-pos", opts_.sgd.sigma_context_pos_max, opts_.sgd.sigma_context_pos_max);
    ops >> Option('d', "sigma-decay", opts_.sgd.sigma_decay, opts_.sgd.sigma_decay);
    ops >> Option('p', "sigma-pc", opts_.sgd.sigma_pc_max, opts_.sgd.sigma_pc_max);
    ops >> Option('q', "sigma-relax-epoch", opts_.sgd.sigma_relax_epoch, opts_.sgd.sigma_relax_epoch);
    ops >> Option(' ', "sigma-relax-steps", opts_.sgd.sigma_relax_steps, opts_.sgd.sigma_relax_steps);
    ops >> Option('s', "context-penalty", opts_.sgd.context_penalty, opts_.sgd.context_penalty);
    ops >> Option('S', "context-penalty-epoch", opts_.sgd.context_penalty_epoch, opts_.sgd.context_penalty_epoch);
    ops >> Option(' ', "context-penalty-steps", opts_.sgd.context_penalty_steps, opts_.sgd.context_penalty_steps);

    ops >> Option('m', "model", opts_.modelfile, opts_.modelfile);
    ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
    ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
    ops >> Option(' ', "pc-init", opts_.pc_init, opts_.pc_init);
    ops >> Option('r', "seed", opts_.sgd.seed, opts_.sgd.seed);
    ops >> Option(' ', "neff-dir", opts_.neff_dir, opts_.neff_dir);
    ops >> Option(' ', "neff-ext", opts_.neff_ext, opts_.neff_ext);
    ops >> Option(' ', "neff-nsamples", opts_.neff_nsamples, opts_.neff_nsamples);
    ops >> Option(' ', "neff-pc", opts_.neff_pc, opts_.neff_pc);

    opts_.Validate();
}

template<class Abc>
void CSSgdApp<Abc>::PrintBanner() const {
    fputs("Optimize CRF weights by stochastic gradient descent.\n", out_);
}

template<class Abc>
void CSSgdApp<Abc>::PrintUsage() const {
    fputs("Usage: cssgd -i <trainset> -j <valset> -o <outfile> [options]\n", out_);
}

template<class Abc>
void CSSgdApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-35s %s\n", "-i, --trainset <file>", "File with training set");
    fprintf(out_, "  %-35s %s\n", "-j, --valset <file>", "File with validation set");
    fprintf(out_, "  %-35s %s\n", "-o, --outfile <file>", 
            "Output file for best CRF on the validation set");
    fprintf(out_, "  %-35s %s\n", "-O, --outfile-tset <file>", 
            "Output file for last CRF on the training set");
    fprintf(out_, "  %-35s %s\n", "-R, --progress <file>",
            "Progress table output (def=stdout)");

    fprintf(out_, "  %-35s %s (def=%zu)\n", "-K, --states [1,inf[",
            "Number of states in CRF to be trained", opts_.nstates);
    fprintf(out_, "  %-35s %s (def=%zu)\n", "-N, --epochs [1,inf[",
            "Maximal number of SGD epochs", opts_.sgd.max_epochs);
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "-t, --toll [0,1]",
            "Log-likelihood change per column for convergence", opts_.sgd.toll);
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "-T, --early-delta [0,inf[",
            "Deviation from the maximal LL on the validation set for early-stopping", opts_.sgd.early_delta);
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "    --min-ll ]-inf,inf[",
            "Minimum log-likelihood on training set in the first epoch", opts_.sgd.min_ll);
    fprintf(out_, "  %-35s %s (def=%zu)\n", "    --min-ll-repeats [0,inf[",
            "Maximal number of repetitions of the first epoch", opts_.sgd.min_ll_repeats);
    fprintf(out_, "\n");

    fprintf(out_, "  %-35s %s (def=%zu)\n", "-E, --eta-mode [1-2]",
            "Mode for updating the learning rate eta", opts_.sgd.eta_mode);
    fprintf(out_, "  %-35s %s\n", "", "1: ALAP3");
    fprintf(out_, "  %-35s %s\n", "", "2: Harmonic function");
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "-e, --eta [0,1]",
            "Initial Learning rate eta in gradient steps", opts_.sgd.eta_init);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-D, --eta-decay [1,inf[",
            "Decay of harmonic function used for updating the learning rate eta", opts_.sgd.eta_decay);
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "    --eta-reinit [0,1]",
            "Learning rate eta for reinitialization", opts_.sgd.eta_reinit);
    fprintf(out_, "  %-35s %s (def=%zu)\n", "    --eta-reinit-num [0,inf[",
            "Number of eta reinitializations", opts_.sgd.eta_reinit_num);
    fprintf(out_, "  %-35s %s (def=%.3g)\n", "    --eta-reinit-delta [0,inf[",
            "Delta of LL change for reinitializing eta", opts_.sgd.eta_reinit_delta);
    fprintf(out_, "  %-35s %s (def=%zu)\n", "-B, --blocks [1,N]",
            "Number of training blocks", opts_.sgd.nblocks);
    fprintf(out_, "  %-35s %s (def=%zu)\n", "-P, --prior [1-3]",
            "Prior of the likelihood function", opts_.prior);
    fprintf(out_, "  %-35s %s\n", "", "1: Lasso prior");
    fprintf(out_, "  %-35s %s\n", "", "2: Gaussian prior");
    fprintf(out_, "  %-35s %s\n", "", "3: Unsymmetric prior");
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-b, --sigma-bias ]0,inf]",
            "Std. deviation sigma in prior of bias weights",
            opts_.sgd.sigma_bias);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-c, --sigma-context ]0,inf]",
            "Std. deviation sigma in prior of context weights",
            opts_.sgd.sigma_context);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-C, --sigma-context-pos ]0,inf]",
            "Std. deviation sigma in unsymmetric prior of positive context weights",
            opts_.sgd.sigma_context_pos_max);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-d, --sigma-decay [0,1]",
            "Exponential decay of sigma of context weights", opts_.sgd.sigma_decay);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-p, --sigma-pc ]0,inf]",
            "Std. deviation sigma in prior of pseudocounts weights",
            opts_.sgd.sigma_pc_max);
    fprintf(out_, "  %-35s %s (def=off)\n", "-q, --sigma-relax-epoch [0,inf[",
            "SGD epoche for beginning to relax the prior");
    fprintf(out_, "  %-35s %s (def=%zu)\n", "    --sigma-relax-steps [0,inf[",
            "Number of epochs for relaxing the prior", opts_.sgd.sigma_relax_steps);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "-s, --context-penalty [0,inf[",
            "Penalty that each context-weights column refers to a density distribution", 
            opts_.sgd.context_penalty);
    fprintf(out_, "  %-35s %s (def=off)\n", "-S, --context-penalty-epoch [0,inf[",
            "SGD epoche for beginning to relax the context penalty");
    fprintf(out_, "  %-35s %s (def=%zu)\n", "    --context-penalty-steps [0,inf[",
            "Number of epochs for relaxing the context penalty", opts_.sgd.context_penalty_steps);
    fprintf(out_, "\n");

    fprintf(out_, "  %-35s %s\n", "-m, --model <file>",
            "Model file with CRF or context library for initialization");
    fprintf(out_, "  %-35s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
           "Weight of central profile column in CRF initialization", opts_.weight_center);
    fprintf(out_, "  %-35s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
            "Parameter for exponential decay of window weights in CRF initialization", 
            opts_.weight_decay);
    fprintf(out_, "  %-35s %s (def=off)\n", "-g, --gauss-init [0,inf[",
            "Turn on gaussian CRF initialization using given sigma");
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "    --pc-init ]0,1]",
            "Constant pseudocount admix in CRF initialization by sampling",
            opts_.pc_init);
    fprintf(out_, "  %-35s %s (def=%u)\n", "-r, --seed [0,inf[",
            "Seed for random number generator", opts_.sgd.seed);
    fprintf(out_, "\n");

    fprintf(out_, "  %-35s %s (def=%s)\n", "    --neff-dir <dir>",
           "Directory with sequences or profiles to be used for calculating the Neff", opts_.neff_dir.c_str());
    fprintf(out_, "  %-35s %s (def=%s)\n", "    --neff-ext <file-ext>",
           "File extension of sequences or profiles to be used for calculating the Neff", opts_.neff_ext.c_str());
    fprintf(out_, "  %-35s %s (def=%zu)\n", "    --neff-nsamples [0,inf[",
           "Number of samples to be used for calculating the Neff", opts_.neff_nsamples);
    fprintf(out_, "  %-35s %s (def=%.2f)\n", "    --neff-pc ]0,1]",
           "Pseudocounts admix for calculating the Neff", opts_.neff_pc);
}

template<class Abc>
int CSSgdApp<Abc>::Run() {
    FILE* fin;
    bool trainProfiles;
    bool valProfiles;

    fin = fopen(opts_.trainfile.c_str(), "r");
    if (!fin)
        throw Exception("Can't read training set from '%s'!", opts_.trainfile.c_str());
    trainProfiles = TrainingProfile<Abc>::IsTrainingProfile(fin);
    fclose(fin);

    if (opts_.valfile.empty()) { 
        valProfiles = trainProfiles;
    } else {
        fin = fopen(opts_.valfile.c_str(), "r");
        if (!fin)
            throw Exception("Can't read validation set from '%s'!", opts_.valfile.c_str());
        valProfiles = TrainingProfile<Abc>::IsTrainingProfile(fin);
        fclose(fin);
    }

    if (trainProfiles && valProfiles)
       return CSSgdRunner<Abc, TrainingProfile<Abc>, TrainingProfile<Abc> >(out_, opts_).Run();
    else if (!trainProfiles && valProfiles)
       return CSSgdRunner<Abc, TrainingSequence<Abc>, TrainingProfile<Abc> >(out_, opts_).Run();
    else if (trainProfiles && !valProfiles)
       return CSSgdRunner<Abc, TrainingProfile<Abc>, TrainingSequence<Abc> >(out_, opts_).Run();
    else
       return CSSgdRunner<Abc, TrainingSequence<Abc>, TrainingSequence<Abc> >(out_, opts_).Run();
}


}  // namespace cs

int main(int argc, char* argv[]) {
    int rv = cs::CSSgdApp<cs::AA>().main(argc, argv, stdout, "cssgd");
    return rv;
}

