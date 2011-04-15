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

struct CSSgdAppOptions {
    CSSgdAppOptions() { Init(); }

    void Init() {
        nstates       = 100;
        blosum_type   = "BLOSUM62";
        pc_init       = 0.5f;
        gauss_init    = 0.0;
        prior         = 2;
        sigma_context = 0.75;
        sigma_decay   = 1.0;
        sigma_bias    = 1.0;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (trainfile.empty()) throw Exception("No training set provided!");
        if (valfile.empty()) throw Exception("No validation set provided!");
        if (crffile.empty()) throw Exception("No output file for CRF provided!");
        if (prior < 1 || prior > 2) throw Exception("Invalid prior!");
    }

    void PrintOptions(FILE* out) const {
        fprintf(out, "  %-15s: %s\n", "--trainset", trainfile.c_str()); 
        fprintf(out, "  %-15s: %s\n", "--valset", valfile.c_str()); 
        fprintf(out, "  %-15s: %s\n", "--outfile", crffile.c_str()); 
        fprintf(out, "  %-15s: %s\n", "--progress", outfile.c_str()); 
        fprintf(out, "  %-15s: %s\n", "--model", modelfile.c_str()); 
        fprintf(out, "  %-15s: %zu\n", "--states", nstates); 
        fprintf(out, "  %-15s: %u\n", "--epochs", sgd.max_epochs); 
        fprintf(out, "  %-15s: %1.5f\n", "--toll", sgd.toll); 
        fprintf(out, "  %-15s: %1.5f\n", "--eta", sgd.eta); 
        fprintf(out, "  %-15s: %u\n", "--seed", sgd.seed); 
        fprintf(out, "  %-15s: %1.5f\n", "--gauss-init", gauss_init); 
        fprintf(out, "  %-15s: %zu\n", "--prior", prior); 
        fprintf(out, "  %-15s: %1.5f\n", "--sigma_context", sigma_context); 
        fprintf(out, "  %-15s: %1.5f\n", "--sigma_bias", sigma_bias); 
        fprintf(out, "  %-15s: %1.5f\n", "--sigma_decay", sigma_decay); 
        fprintf(out, "  %-15s: %1.5f\n", "--pc-init", pc_init); 
        fprintf(out, "  %-15s: %zu\n", "--blocks", sgd.nblocks); 
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
    // Sigma of gaussian prior for context weights
    double sigma_context;
    // Parameter governing exponential decay of context sigma
    double sigma_decay;
    // Sigma of gaussian prior for bias weights
    double sigma_bias;
    // Wrapper for SGD parameters
    SgdParams sgd;
};  // CSSgdAppOptions



template<class Abc, class TrainingPair>
struct CSSgdRunner {
    typedef vector<TrainingPair> TrainingSet;

    CSSgdRunner(FILE* const out, const CSSgdAppOptions& opts) : 
        out_(out),
        opts_(opts) 
    {}


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
        fprintf(out_, "%zu records read\n", trainset_.size());

        if (!opts_.valfile.empty()) {
            fprintf(out_, "Reading validation set from %s ...\n",
                    GetBasename(opts_.valfile).c_str());
            fin = fopen(opts_.valfile.c_str(), "r");
            if (!fin)
                throw Exception("Can't read validation set from '%s'!",
                                opts_.valfile.c_str());
            ReadAll(fin, valset_);
            fclose(fin);
            fprintf(out_, "%zu records read\n", valset_.size());
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

            TransformToLin(lib);
            LibraryBasedCrfInit<Abc> init(lib, *sm_);
            size_t nstates = opts_.nstates == 0 ? lib.size() : opts_.nstates;
            crf_.reset(new Crf<Abc>(nstates, lib.wlen(), init));

        } else if (opts_.gauss_init != 0.0) {  // init by sampling weights from gaussian
            fputs("Initializing CRF by sampling weights from gaussian ...\n", out_);
            GaussianCrfInit<Abc> init(opts_.gauss_init, *sm_, opts_.sgd.seed);
            crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));

        } else {  // init by sampling from training set
            fputs("Initializing CRF by sampling from training set ...\n", out_);
            MatrixPseudocounts<Abc> pc(*sm_);
            ConstantAdmix admix(opts_.pc_init);
            SamplingCrfInit<Abc, TrainingPair> init(trainset_, pc, admix, *sm_, opts_.sgd.seed);
            crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));
        }
    }

    // Runs stochastic gradient descent.
    int Run() {
        fputs("Running stochastic gradient descent with following parameters:\n", out_);
        opts_.PrintOptions(out_);
        fputs("\n", out_);

        InitSubstitutionMatrix();
        ReadTrainingData();
        InitCrf();

        FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");
        fputs("\n", fout);

        CrfFunc<Abc, TrainingPair> val_func(valset_, *sm_);
       
        DerivCrfFuncPrior<Abc>* prior;
        if (opts_.prior == 1) {
            prior = new GaussianDerivCrfFuncPrior<Abc>(
                opts_.sigma_context, 
                opts_.sigma_decay, 
                opts_.sigma_bias);
        } else {
            prior = new LassoDerivCrfFuncPrior<Abc>(
                opts_.sigma_context, 
                opts_.sigma_decay, 
                opts_.sigma_bias);
        }
        DerivCrfFunc<Abc, TrainingPair> train_func(trainset_, *sm_, *prior);
        SgdOptimizer<Abc, TrainingPair> sgd(train_func, val_func, opts_.sgd);
        sgd.Optimize(*crf_, fout);
        delete prior;

        if (!opts_.outfile.empty()) fclose(fout);

        fout = fopen(opts_.crffile.c_str(), "w");
        if (!fout) throw Exception("Can't write to file '%s'!", opts_.crffile.c_str());
        crf_->Write(fout);
        fclose(fout);
        fprintf(out_, "\nWrote best CRF on validation set to %s\n", opts_.crffile.c_str());

        return 0;
    }

    FILE* const out_;
    const CSSgdAppOptions opts_;
    TrainingSet trainset_;
    TrainingSet valset_;
    scoped_ptr<Crf<Abc> > crf_;
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;

};  //CSSgdRunner



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
int CSSgdApp<Abc>::Run() {
    FILE* fin;
    fin = fopen(opts_.trainfile.c_str(), "r");
    if (!fin)
        throw Exception("Can't read training set from '%s'!", opts_.trainfile.c_str());
    if (StreamStartsWith(fin, "TrainingSequence")) {
        return CSSgdRunner<Abc, TrainingSequence<Abc> >(out_, opts_).Run();
    } else {
        rewind(fin);
        if (StreamStartsWith(fin, "TrainingProfile"))
            return CSSgdRunner<Abc, TrainingProfile<Abc> >(out_, opts_).Run(); 
        else
            throw Exception("Training set does not start with class id 'TrainingSequence' or 'TrainingProfile'");
    }
}

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
    ops >> Option('e', "eta", opts_.sgd.eta, opts_.sgd.eta);
    ops >> Option('r', "seed", opts_.sgd.seed, opts_.sgd.seed);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
    ops >> Option('P', "prior", opts_.prior, opts_.prior);
    ops >> Option('c', "sigma-context", opts_.sigma_context, opts_.sigma_context);
    ops >> Option('b', "sigma-bias", opts_.sigma_bias, opts_.sigma_bias);
    ops >> Option('d', "sigma-decay", opts_.sigma_decay, opts_.sigma_decay);
    ops >> Option('p', "pc-init", opts_.pc_init, opts_.pc_init);
    ops >> Option('B', "blocks", opts_.sgd.nblocks, opts_.sgd.nblocks);

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
    fprintf(out_, "  %-30s %s (def=%6.2g)\n", "-t, --toll [0,1]",
            "Log-likelihood change per column for convergence", opts_.sgd.toll);
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
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-c, --sigma-context [0,inf]",
            "Std. deviation sigma in prior of context weights",
            opts_.sigma_context);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-b, --sigma-bias [0,inf]",
            "Std. deviation sigma in prior of bias weights",
            opts_.sigma_bias);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-d, --sigma-decay [0,1]",
            "Exponential decay of sigma for context weights", opts_.sigma_decay);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-p, --pc-init [0,1]",
            "Constant pseudocount admix in CRF initialization by sampling",
            opts_.pc_init);
}

}  // namespace cs

int main(int argc, char* argv[]) {
    int rv = cs::CSSgdApp<cs::AA>().main(argc, argv, stdout, "cssgd");
    return rv;
}

