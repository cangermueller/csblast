// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-31		

#include "cs.h"
#include "application.h"
#include "blosum_matrix.h"
#include "getopt_pp.h"
#include "lr_sgd-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"

using namespace GetOpt;
using std::string;
using std::vector;


namespace cs {


// CSSgdAppOptions


// Stores global application options.
struct CSSgdAppOptions {
    CSSgdAppOptions() { Init(); }

    void Init() {
		wlen		  = 13;
		klen		  = 4;
        blosum_type   = "BLOSUM62";
        gauss_init    = 0.1;
        sigma_bias    = 1.0;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (trainfile.empty()) throw Exception("No training set provided!");
        if (valfile.empty()) throw Exception("No validation set provided!");
        if (lr_file.empty()) throw Exception("No output file for lr params provided!");
    }

	void PrintOptions(FILE* out) const {
		fprintf(out, "  %-15s: %2zu\n", "--wlen", wlen); 
		fprintf(out, "  %-15s: %2zu\n", "--klen", klen); 
		fprintf(out, "  %-15s: %s\n", "--trainfile", trainfile.c_str()); 
		fprintf(out, "  %-15s: %s\n", "--valfile", valfile.c_str()); 
		if (!modelfile.empty())
			fprintf(out, "  %-15s: %s\n", "modelfile", modelfile.c_str()); 
		if (!outfile.empty())
			fprintf(out, "  %-15s: %s\n", "outfile", outfile.c_str()); 
		fprintf(out, "  %-15s: %s\n", "blosum type", blosum_type.c_str()); 
		if (modelfile.empty()) 
			fprintf(out, "  %-15s: %5.2f\n", "sigma gauss", gauss_init); 
		fprintf(out, "  %-15s: %5.2f\n", "sigma bias", sigma_bias); 
	}

	// Variables
    string trainfile;	// Input file with training set.
    string valfile;		// Input file with validation.
	string modelfile;	// Input file for lr params initialization.
    string lr_file;		// Output file for lr parameters
    string outfile;		// Output file for progress table.
	size_t wlen;		// Window length
	size_t klen;		// Kmer length
    string blosum_type;	// BLOSUM matrix for pseudocount generation.
    double gauss_init;	// Sigma for lr init from gaussian
    double sigma_bias;	// Sigma of gaussian prior for bias weights
    LrSgdParams sgd;	// Wrapper for SGD parameters
};  // CSSgdAppOptions


// CSSgdApp


// Main application.
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
    // Initializes lr params by sampling the gaussian distribution or using a file for jumpstart.
    void InitLrParams();
    // Initializes substitution matrix (specialized by alphabet type).
    void InitSubstitutionMatrix();


	// Variables
    CSSgdAppOptions opts_;						// global application options
    TrainingSet trainset_;						// training set
    TrainingSet valset_;						// validation set
	scoped_ptr<LrCfg<Abc> > cfg_;				// window configs
    scoped_ptr<LrParams<Abc> > lr_params_;		// lr params being optimized
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;	// substitution matrix
};  // CSSgdApp


template<class Abc>
void CSSgdApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "trainset", opts_.trainfile, opts_.trainfile);
    ops >> Option('j', "valset", opts_.valfile, opts_.valfile);
    ops >> Option('o', "outfile", opts_.lr_file, opts_.lr_file);
    ops >> Option('O', "progress", opts_.outfile, opts_.outfile);
    ops >> Option('m', "model", opts_.modelfile, opts_.modelfile);
    ops >> Option('w', "wlen", opts_.wlen, opts_.wlen);
    ops >> Option('k', "klen", opts_.klen, opts_.klen);
    ops >> Option('N', "epochs", opts_.sgd.max_epochs, opts_.sgd.max_epochs);
    ops >> Option('t', "toll", opts_.sgd.toll, opts_.sgd.toll);
    ops >> Option('e', "eta", opts_.sgd.eta, opts_.sgd.eta);
    ops >> Option('r', "seed", opts_.sgd.seed, opts_.sgd.seed);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
    ops >> Option('b', "sigma-bias", opts_.sigma_bias, opts_.sigma_bias);
    ops >> Option('B', "blocks", opts_.sgd.nblocks, opts_.sgd.nblocks);

    opts_.Validate();
}

template<class Abc>
void CSSgdApp<Abc>::PrintBanner() const {
    fputs("Optimize lr params by stochastic gradient descent.\n", out_);
}

template<class Abc>
void CSSgdApp<Abc>::PrintUsage() const {
    fputs("Usage: cssgd -i <trainset> -j <valset> -o <outfile> [options]\n", out_);
}

template<class Abc>
void CSSgdApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n", "-i, --trainset <file>", "File with training set");
    fprintf(out_, "  %-30s %s\n", "-j, --valset <file>", "File with validation set");
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "lr params output file");
    fprintf(out_, "  %-30s %s\n", "-O, --progress <file>",
            "Progress table output (def=stdout)");
    fprintf(out_, "  %-30s %s\n", "-m, --model <file>",
            "Model file with lr params for initialization");
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-w, --wlen [1,inf[",
            "Window length of training sequences", opts_.wlen);
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-k, --klen [1,inf[",
            "Length of kmers within a window", opts_.klen);
    fprintf(out_, "  %-30s %s (def=%d)\n", "-N, --epochs [1,inf[",
            "Maximal number of SGD epochs", opts_.sgd.max_epochs);
    fprintf(out_, "  %-30s %s (def=%6.2g)\n", "-t, --toll [0,1]",
            "Log-likelihood change per column for convergence", opts_.sgd.toll);
    fprintf(out_, "  %-30s %s (def=%.3f)\n", "-e, --eta [0,1]",
            "Learning rate eta in gradient steps", opts_.sgd.eta);
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-B, --blocks [1,N]",
            "Number of training blocks", opts_.sgd.nblocks);
    fprintf(out_, "  %-30s %s (def=off)\n", "-g, --gauss-init [0,inf[",
            "Turn on gaussian lr params initialization using given sigma");
    fprintf(out_, "  %-30s %s (def=%u)\n", "-r, --seed [0,inf[",
            "Seed for random number generator", opts_.sgd.seed);
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-b, --sigma-bias [0,inf]",
            "Std. deviation sigma in gaussian prior of bias weights",
            opts_.sigma_bias);
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
void CSSgdApp<Abc>::InitLrParams() {
	cfg_.reset(new LrCfg<Abc>(opts_.wlen, opts_.klen));
    if (!opts_.modelfile.empty() && GetFileExt(opts_.modelfile) == "lr") {
		// init by modelfile
        fprintf(out_, "Reading lr params from %s ...\n", GetBasename(opts_.modelfile).c_str());

        FILE* fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        lr_params_.reset(new LrParams<Abc>(fin));
        fclose(fin);

    } else {  
		// init by sampling from gaussian
        fputs("Initializing lr params by sampling from gaussian ...\n", out_);
        GaussianLrParamsInit<Abc> init(opts_.gauss_init, opts_.sgd.seed);
        lr_params_.reset(new LrParams<Abc>(*cfg_, init));
    }
}

template<class Abc>
int CSSgdApp<Abc>::Run() {
	fputs("Running stochastic gradient descent with following parameters:\n", out_);
	opts_.PrintOptions(out_);
    fputs("\n", out_);

    InitSubstitutionMatrix();
    ReadTrainingData();
    InitLrParams();

    FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");
    fputs("\n", fout);

    // Run stochastic gradient descent
    LrFunc<Abc, TrainingSequence<Abc> > val_func(*cfg_, valset_, *sm_);
    DerivLrFunc<Abc, TrainingSequence<Abc> > train_func(*cfg_, trainset_, *sm_,
                                                         opts_.sigma_bias);
    LrSgdOptimizer<Abc, TrainingSequence<Abc> > sgd(train_func, val_func, opts_.sgd);
    sgd.Optimize(*lr_params_, fout);

    if (!opts_.outfile.empty()) fclose(fout);

    fout = fopen(opts_.lr_file.c_str(), "w");
    if (!fout) throw Exception("Can't write to file '%s'!", opts_.lr_file.c_str());
    lr_params_->Write(fout);
    fclose(fout);
    fprintf(out_, "\nWrote best lr params on validation set to %s\n", opts_.lr_file.c_str());

    return 0;
}


}  // namespace cs





int main(int argc, char* argv[]) {
    int rv = cs::CSSgdApp<cs::AA>().main(argc, argv, stdout, "cssgd");
    return rv;
}

