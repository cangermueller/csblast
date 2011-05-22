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
        sigma_bias     = 10.0;
        sigma_context  = 10.0;
        sigma_decay    = 1.0;
        weight_center  = 1.6;
        weight_decay   = 0.85;
        neff_ext       = "seq";
        neff_nsamples  = 100;
        neff_pc        = 1.0;
        save           = false;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (trainfile.empty()) throw Exception("No training set provided!");
        if (valfile.empty()) throw Exception("No validation set provided!");
        if (crffile.empty()) throw Exception("No output file for CRF provided!");
        if (prior < 1 || prior > 2) throw Exception("Invalid prior!");
        if (pc_init <= 0 || pc_init > 1.0) throw Exception("Pseudocounts admix for initialization invalid!");
        if (neff_pc <= 0 || neff_pc > 1.0) throw Exception("Pseudocounts admix for computing the Neff invalid!");
    }

    void PrintOptions(FILE* out) const {
        fprintf(out, "  %-20s: %s\n", "--trainset", trainfile.c_str()); 
        fprintf(out, "  %-20s: %s\n", "--valset", valfile.c_str()); 
        fprintf(out, "  %-20s: %s\n", "--outfile", crffile.c_str()); 
        fprintf(out, "  %-20s: %s\n", "--progress", outfile.c_str()); 
        fprintf(out, "  %-20s: %s\n", "--model", modelfile.c_str()); 
        fprintf(out, "  %-20s: %zu\n", "--states", nstates); 
        fprintf(out, "  %-20s: %u\n", "--epochs", sgd.max_epochs); 
        fprintf(out, "  %-20s: %1.5f\n", "--toll", sgd.toll); 
        fprintf(out, "  %-20s: %1.5f\n", "--early-delta", sgd.early_delta); 
        fprintf(out, "  %-20s: %1.5f\n", "--eta", sgd.eta); 
        fprintf(out, "  %-20s: %u\n", "--seed", sgd.seed); 
        fprintf(out, "  %-20s: %1.5f\n", "--gauss-init", gauss_init); 
        fprintf(out, "  %-20s: %zu\n", "--prior", prior); 
        fprintf(out, "  %-20s: %1.5f\n", "--sigma_bias", sigma_bias); 
        fprintf(out, "  %-20s: %1.5f\n", "--sigma_context", sigma_context); 
        fprintf(out, "  %-20s: %1.5f\n", "--sigma_decay", sigma_decay); 
        fprintf(out, "  %-20s: %1.5f\n", "--sigma_pc", sgd.sigma_pc_max); 
        fprintf(out, "  %-20s: %zu\n", "--sigma_pc_epoch", sgd.sigma_pc_epoch); 
        fprintf(out, "  %-20s: %1.5f\n", "--sigma_pc_delta", sgd.sigma_pc_delta); 
        fprintf(out, "  %-20s: %1.5f\n", "--pc-init", pc_init); 
        fprintf(out, "  %-20s: %zu\n", "--blocks", sgd.nblocks); 
        fprintf(out, "  %-20s: %.2f\n", "--weight_center", weight_center); 
        fprintf(out, "  %-20s: %.2f\n", "--weight_decay", weight_decay); 
        fprintf(out, "  %-20s: %s\n", "--neff-dir", neff_dir.c_str()); 
        fprintf(out, "  %-20s: %s\n", "--neff-ext", neff_ext.c_str()); 
        fprintf(out, "  %-20s: %zu\n", "--neff-nsamples", neff_nsamples); 
        fprintf(out, "  %-20s: %.2f\n", "--neff-pc", neff_pc); 
        fprintf(out, "  %-20s: %d\n", "--save", save); 
    }


    // Input file with training set.
    string trainfile;
    // Input file with validation.
    string valfile;
    // Output file for CRF
    string crffile;
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
    // Sigma of prior for context weights
    double sigma_context;
    // Parameter governing exponential decay of context sigma
    double sigma_decay;
    // Sigma of prior for bias weights
    double sigma_bias;
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
    // Save currently best CRF on-line
    bool save;

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
            prior.reset(new GaussianDerivCrfFuncPrior<Abc>(
                opts_.sigma_bias,
                opts_.sigma_context, 
                opts_.sigma_decay, 
                opts_.sgd.sigma_pc_max));
        else
            prior.reset(new LassoDerivCrfFuncPrior<Abc>(
                opts_.sigma_bias,
                opts_.sigma_context, 
                opts_.sigma_decay, 
                opts_.sgd.sigma_pc_max));

        CrfFunc<Abc, TrainingPairV> val_func(valset_, *sm_);
        DerivCrfFunc<Abc, TrainingPairT> train_func(trainset_, *sm_, *prior);
        SgdOptimizer<Abc, TrainingPairT, TrainingPairV> sgd(train_func, val_func, 
                opts_.sgd, neff_samples_, opts_.neff_pc, opts_.save ? opts_.crffile : "");
        sgd.Optimize(*crf_, fout);

        if (!opts_.outfile.empty()) fclose(fout);

        if (!opts_.save) {
            fout = fopen(opts_.crffile.c_str(), "w");
            if (!fout) throw Exception("Can't write to file '%s'!", opts_.crffile.c_str());
            crf_->Write(fout);
            fclose(fout);
        }
        fprintf(out_, "\nWrote CRF to %s!\n", opts_.crffile.c_str());

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
    ops >> Option('o', "outfile", opts_.crffile, opts_.crffile);
    ops >> Option('O', "progress", opts_.outfile, opts_.outfile);
    ops >> Option('m', "model", opts_.modelfile, opts_.modelfile);
    ops >> Option('K', "states", opts_.nstates, opts_.nstates);
    ops >> Option('N', "epochs", opts_.sgd.max_epochs, opts_.sgd.max_epochs);
    ops >> Option('t', "toll", opts_.sgd.toll, opts_.sgd.toll);
    ops >> Option('T', "early-delta", opts_.sgd.early_delta, opts_.sgd.early_delta);
    ops >> Option('e', "eta", opts_.sgd.eta, opts_.sgd.eta);
    ops >> Option('r', "seed", opts_.sgd.seed, opts_.sgd.seed);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
    ops >> Option('P', "prior", opts_.prior, opts_.prior);
    ops >> Option('p', "sigma-pc", opts_.sgd.sigma_pc_max, opts_.sgd.sigma_pc_max);
    ops >> Option('q', "sigma-pc-epoch", opts_.sgd.sigma_pc_epoch, opts_.sgd.sigma_pc_epoch);
    ops >> Option('Q', "sigma-pc-delta", opts_.sgd.sigma_pc_delta, opts_.sgd.sigma_pc_delta);
    ops >> Option('c', "sigma-context", opts_.sigma_context, opts_.sigma_context);
    ops >> Option('d', "sigma-decay", opts_.sigma_decay, opts_.sigma_decay);
    ops >> Option('b', "sigma-bias", opts_.sigma_bias, opts_.sigma_bias);
    ops >> Option('B', "blocks", opts_.sgd.nblocks, opts_.sgd.nblocks);
    ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
    ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
    ops >> Option(' ', "neff-dir", opts_.neff_dir, opts_.neff_dir);
    ops >> Option(' ', "neff-ext", opts_.neff_ext, opts_.neff_ext);
    ops >> Option(' ', "neff-nsamples", opts_.neff_nsamples, opts_.neff_nsamples);
    ops >> Option(' ', "neff-pc", opts_.neff_pc, opts_.neff_pc);
    ops >> Option(' ', "pc-init", opts_.pc_init, opts_.pc_init);
    ops >> OptionPresent(' ', "save", opts_.save);

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
    fprintf(out_, "  %-30s %s\n", "-i, --trainset <file>", "File with training set");
    fprintf(out_, "  %-30s %s\n", "-j, --valset <file>", "File with validation set");
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "CRF output file");
    fprintf(out_, "  %-30s %s\n", "-O, --progress <file>",
            "Progress table output (def=stdout)");
    fprintf(out_, "  %-30s %s\n", "-m, --model <file>",
            "Model file with CRF or context library for initialization");
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-K, --states [1,inf[",
            "Number of states in CRF to be trained", opts_.nstates);
    fprintf(out_, "  %-30s %s (def=%d)\n", "-N, --epochs [1,inf[",
            "Maximal number of SGD epochs", opts_.sgd.max_epochs);
    fprintf(out_, "  %-30s %s (def=%.2g)\n", "-t, --toll [0,1]",
            "Log-likelihood change per column for convergence", opts_.sgd.toll);
    fprintf(out_, "  %-30s %s (def=%.2g)\n", "-T, --early-delta [0,inf[",
            "Deviation from the maximal LL on the validation set for early-stopping", opts_.sgd.early_delta);
    fprintf(out_, "  %-30s %s (def=%.3f)\n", "-e, --eta [0,1]",
            "Learning rate eta in gradient steps", opts_.sgd.eta);
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-B, --blocks [1,N]",
            "Number of training blocks", opts_.sgd.nblocks);
    fprintf(out_, "  %-30s %s (def=off)\n", "-g, --gauss-init [0,inf[",
            "Turn on gaussian CRF initialization using given sigma");
    fprintf(out_, "  %-30s %s (def=%u)\n", "-r, --seed [0,inf[",
            "Seed for random number generator", opts_.sgd.seed);
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-P, --prior [1-2]",
            "Prior of the likelihood function", opts_.prior);
    fprintf(out_, "  %-30s %s\n", "", "1: gaussian prior");
    fprintf(out_, "  %-30s %s\n", "", "2: lasso prior");
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-p, --sigma-pc ]0,inf]",
            "Std. deviation sigma in prior of pseudocounts weights",
            opts_.sgd.sigma_pc_max);
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-q, --sigma-pc-epoch [0,inf[",
            "SGD epoche for activating prior of pseudocounts weights", opts_.sgd.sigma_pc_epoch);
    fprintf(out_, "  %-30s %s (def=%.3f)\n", "-r, --sigma-pc-delta [0,inf[",
            "Gradient for activating prior of pseudocounts weights", opts_.sgd.sigma_pc_delta);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-c, --sigma-context ]0,inf]",
            "Std. deviation sigma in prior of context weights",
            opts_.sigma_context);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-d, --sigma-decay [0,1]",
            "Exponential decay of sigma for context weights", opts_.sigma_decay);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-b, --sigma-bias ]0,inf]",
            "Std. deviation sigma in prior of bias weights",
            opts_.sigma_bias);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --pc-init ]0,1]",
            "Constant pseudocount admix in CRF initialization by sampling",
            opts_.pc_init);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
           "Weight of central profile column in CRF initialization", opts_.weight_center);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
            "Parameter for exponential decay of window weights in CRF initialization", 
            opts_.weight_decay);
    fprintf(out_, "  %-30s %s (def=%s)\n", "    --neff-dir <dir>",
           "Directory with sequences or profiles to be used for calculating the Neff", opts_.neff_dir.c_str());
    fprintf(out_, "  %-30s %s (def=%s)\n", "    --neff-ext <file-ext>",
           "File extension of sequences or profiles to be used for calculating the Neff", opts_.neff_ext.c_str());
    fprintf(out_, "  %-30s %s (def=%zu)\n", "    --neff-nsamples [0,inf[",
           "Number of samples to be used for calculating the Neff", opts_.neff_nsamples);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --neff-pc ]0,1]",
           "Pseudocounts admix for calculating the Neff", opts_.neff_pc);
    fprintf(out_, "  %-30s %s (def=off)\n", "    --save",
           "Save currently best CRF on-line");
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

