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
        sigma_context = 0.2;
        sigma_decay   = 0.9;
        sigma_bias    = 1.0;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (trainfile.empty()) throw Exception("No training set provided!");
        if (valfile.empty()) throw Exception("No validation set provided!");
        if (crffile.empty()) throw Exception("No output file for CRF provided!");
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
    // Sigma of gaussian prior for context weights
    double sigma_context;
    // Parameter governing exponential decay of context sigma
    double sigma_decay;
    // Sigma of gaussian prior for bias weights
    double sigma_bias;
    // Wrapper for SGD parameters
    SgdParams sgd;
};  // CSSgdAppOptions


template<class Abc>
class CSSgdApp : public Application {
  protected:
    typedef vector<TrainingSequence<Abc> > TrainingSet;

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
    // Reads training and validation set.
    void ReadTrainingData();
    // Initializes CRF from context library, training data, or jumpstart file.
    void InitCrf();
    // Initializes substitution matrix (specialized by alphabet type).
    void InitSubstitutionMatrix();

    CSSgdAppOptions opts_;
    TrainingSet trainset_;
    TrainingSet valset_;
    scoped_ptr<Crf<Abc> > crf_;
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;
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
    ops >> Option('e', "eta", opts_.sgd.eta, opts_.sgd.eta);
    ops >> Option('r', "seed", opts_.sgd.seed, opts_.sgd.seed);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
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
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-c, --sigma-context [0,inf]",
            "Std. deviation sigma in gaussian prior of context weights",
            opts_.sigma_context);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-b, --sigma-bias [0,inf]",
            "Std. deviation sigma in gaussian prior of bias weights",
            opts_.sigma_bias);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-d, --sigma-decay [0,1]",
            "Exponential decay of sigma for context weights", opts_.sigma_decay);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-p, --pc-init [0,1]",
            "Constant pseudocount admix in CRF initialization by sampling",
            opts_.pc_init);
}

template<class Abc>
void CSSgdApp<Abc>::InitSubstitutionMatrix() {
    BlosumType type = BlosumTypeFromString(opts_.blosum_type);
    sm_.reset(new BlosumMatrix(type));
}

template<class Abc>
void CSSgdApp<Abc>::ReadTrainingData() {
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

template<class Abc>
void CSSgdApp<Abc>::InitCrf() {
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
        SamplingCrfInit<Abc, TrainingSequence<Abc> > init(trainset_, pc, admix,
                                                          *sm_, opts_.sgd.seed);
        crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));
    }
}

template<class Abc>
int CSSgdApp<Abc>::Run() {
    InitSubstitutionMatrix();
    ReadTrainingData();
    InitCrf();

    FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");
    fputs("\n", fout);

    // Run stochastic gradient descent
    CrfFunc<Abc, TrainingSequence<Abc> > val_func(valset_, *sm_);
    DerivCrfFunc<Abc, TrainingSequence<Abc> > train_func(trainset_, *sm_,
                                                         opts_.sigma_context,
                                                         opts_.sigma_decay,
                                                         opts_.sigma_bias);
    SgdOptimizer<Abc, TrainingSequence<Abc> > sgd(train_func, val_func, opts_.sgd);
    sgd.Optimize(*crf_, fout);

    if (!opts_.outfile.empty()) fclose(fout);

    fout = fopen(opts_.crffile.c_str(), "w");
    if (!fout) throw Exception("Can't write to file '%s'!", opts_.crffile.c_str());
    crf_->Write(fout);
    fclose(fout);
    fprintf(out_, "\nWrote best CRF on validation set to %s\n", opts_.crffile.c_str());

    return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
    int rv = cs::CSSgdApp<cs::AA>().main(argc, argv, stdout, "cssgd");
    return rv;
}

