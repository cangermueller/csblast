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

#include "cs.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf-inl.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "hmc.h"
#include "matrix_pseudocounts-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct CSHmcAppOptions {
  CSHmcAppOptions() { Init(); }

  void Init() {
    nstates        = 100;
    nsteps         = 100;
    blosum_type    = "BLOSUM62";
    pc_init        = 0.5f;
    gauss_init     = 0.0;
    theta          = 1.0;
    prior          = 1;
    sigma_context  = 0.2;
    sigma_decay    = 0.9;
    sigma_bias     = 1.0;
    no_sort        = false;
    hmc.sgd_epochs = 0;
    sgd.eta        = 0.001;
    weight_center  = 1.6;
    weight_decay   = 0.85;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (hmc.lfsteps % (2 * hmc.nblocks) != 0)
      throw Exception("Leapfrog steps must be multiple of 2x number of blocks!");
    if (trainfile.empty()) throw Exception("No training set provided!");
    if (valfile.empty()) throw Exception("No validation set provided!");
    if (crffile.empty()) throw Exception("No output file for CRF provided!");
  }

  // Input file with training set.
  string trainfile;
  // Input file with validation.
  string valfile;
  // Output file for CRF.
  string crffile;
  // Output file for progress table.
  string outfile;
  // Crf input file with context library or CRF for initialization.
  string modelfile;
  // The number of states in the HMM to train.
  size_t nstates;
  // Number of HMC sampling steps
  size_t nsteps;
  // BLOSUM matrix for pseudocount generation.
  string blosum_type;
  // Constant admix for CRF initialization
  double pc_init;
  // Sigma for CRF init from gaussian
  double gauss_init;
  // Probability for trying a replica exchange
  double theta;
  // Prior to be used for calculating the likelihood
  size_t prior;
  // Sigma of gaussian prior for context weights
  double sigma_context;
  // Parameter governing exponential decay of context sigma
  double sigma_decay;
  // Sigma of gaussian prior for bias weights
  double sigma_bias;
  // Use sorted parallel tempering instead of single exchange parallel tempering
  bool no_sort;
  // Wrapper for HMC parameters
  HmcParams hmc;
  // Parameter wrapper for SGD in basin hopping
  SgdParams sgd;
  // Weight of central column in multinomial emission
  double weight_center;
  // Exponential decay of window weights
  double weight_decay;
};  // CSHmcAppOptions


template<class Abc>
class CSHmcApp : public Application {
 protected:
  typedef std::vector<TrainingSequence<Abc> > TrainingSet;
  typedef CrfFunc<Abc, TrainingSequence<Abc> > Likelihood;
  typedef DerivCrfFunc<Abc, TrainingSequence<Abc> > Gradient;
  typedef LeapfrogProposal<Abc, TrainingSequence<Abc> > Leapfrog;
  typedef BasinHoppingProposal<Abc, TrainingSequence<Abc> > BasinHopping;
  typedef ParallelTempering<Abc, TrainingSequence<Abc> > PT;
#ifdef PARALLEL
  typedef SingleExchParallelTempering<Abc, TrainingSequence<Abc> > SingleExchPT;
  typedef SortedParallelTempering<Abc, TrainingSequence<Abc> > SortedPT;
#endif

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

  CSHmcAppOptions opts_;
  TrainingSet trainset_;
  TrainingSet valset_;
  scoped_ptr<Crf<Abc> > crf_;
  scoped_ptr<SubstitutionMatrix<Abc> > sm_;
  scoped_ptr<Leapfrog> propose_;
  scoped_ptr<PT> partemp_;
};  // CSHmcApp


template<class Abc>
void CSHmcApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "trainset", opts_.trainfile, opts_.trainfile);
  ops >> Option('j', "valset", opts_.valfile, opts_.valfile);
  ops >> Option('o', "outfile", opts_.crffile, opts_.crffile);
  ops >> Option('O', "progress", opts_.outfile, opts_.outfile);
  ops >> Option('m', "model", opts_.modelfile, opts_.modelfile);
  ops >> Option('K', "states", opts_.nstates, opts_.nstates);
  ops >> Option('N', "steps", opts_.nsteps, opts_.nsteps);
  ops >> Option('L', "leapfrog", opts_.hmc.lfsteps, opts_.hmc.lfsteps);
  ops >> Option('T', "temp", opts_.hmc.temp, opts_.hmc.temp);
  ops >> Option('e', "epsilon", opts_.hmc.epsilon, opts_.hmc.epsilon);
  ops >> Option('r', "seed", opts_.hmc.seed, opts_.hmc.seed);
  ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
  ops >> Option('E', "epochs", opts_.hmc.sgd_epochs, opts_.hmc.sgd_epochs);
  ops >> Option('P', "prior", opts_.prior, opts_.prior);
  ops >> Option('c', "sigma-context", opts_.sigma_context, opts_.sigma_context);
  ops >> Option('b', "sigma-bias", opts_.sigma_bias, opts_.sigma_bias);
  ops >> Option('d', "sigma-decay", opts_.sigma_decay, opts_.sigma_decay);
  ops >> Option('p', "pc-init", opts_.pc_init, opts_.pc_init);
  ops >> Option('t', "theta", opts_.theta, opts_.theta);
  ops >> Option('a', "alpha", opts_.hmc.alpha, opts_.hmc.alpha);
  ops >> Option('B', "blocks", opts_.hmc.nblocks, opts_.hmc.nblocks);
  ops >> Option('S', "sgd-blocks", opts_.sgd.nblocks, opts_.sgd.nblocks);
  ops >> Option(' ', "eps-up", opts_.hmc.epsilon_up, opts_.hmc.epsilon_up);
  ops >> Option(' ', "eps-down", opts_.hmc.epsilon_down, opts_.hmc.epsilon_down);
  ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
  ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
  ops >> OptionPresent(' ', "no-sort", opts_.no_sort);

  opts_.Validate();
}

template<class Abc>
void CSHmcApp<Abc>::PrintBanner() const {
  fputs("Sample CRF parameters with hybrid Monte Carlo.\n", out_);
}

template<class Abc>
void CSHmcApp<Abc>::PrintUsage() const {
  fputs("Usage: cshmc -i <trainset> -j <valset> -o <outfile> [options]\n", out_);
}

template<class Abc>
void CSHmcApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --trainset <file>", "File with training set");
  fprintf(out_, "  %-30s %s\n", "-j, --valset <file>", "File with validation set");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "CRF output file");
  fprintf(out_, "  %-30s %s\n", "-O, --progress <file>",
          "Progress table output (def=stdout)");
  fprintf(out_, "  %-30s %s\n", "-m, --model <file>",
          "Model file with CRF or context library for initialization");
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-K, --states [1,inf[",
          "Number of states in CRF to be trained", opts_.nstates);
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-N, --steps [1,inf[",
          "Number of HMC sampling steps", opts_.nsteps);
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-L, --leapfrog [1,inf[",
          "Number of leapfrog steps", opts_.hmc.lfsteps);
  fprintf(out_, "  %-30s %s\n", "-E, --epochs [1,inf[",
          "Enable basin hopping by specifying number of SGD epochs (def=off)");
  fprintf(out_, "  %-30s %s (def=%6.2g)\n", "-e, --epsilon [0,1]",
          "Start value for time step epsilon in leapfrog", opts_.hmc.epsilon);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "-T, --temp [0,inf[",
          "Temperature level for sampling", opts_.hmc.temp);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "-a, --alpha [1,inf]",
          "Geometric increase of temperatures in adjacent replicas",
          opts_.hmc.alpha);
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-B, --blocks [1,N]",
          "Number of training blocks in online training", opts_.hmc.nblocks);

  fprintf(out_, "  %-30s %s (def=off)\n", "-g, --gauss-init [0,inf[",
          "Turn on gaussian CRF initialization using given sigma");
#ifdef PARALLEL
  fprintf(out_, "  %-30s %s (def=%.3f)\n", "-t, --theta [0,1]",
          "Probability of replica exchange", opts_.theta);
#endif
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-S, --sgd-blocks [1,N]",
          "Number of training blocks in SGD", opts_.sgd.nblocks);
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-P, --prior [1-2]",
          "Prior of the likelihood function", opts_.prior);
  fprintf(out_, "  %-30s %s\n", "", "1: gaussian prior");
  fprintf(out_, "  %-30s %s\n", "", "2: lasso prior");
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
  fprintf(out_, "  %-30s %s (def=%u)\n", "-r, --seed [0,inf[",
          "Seed for random number generator", opts_.hmc.seed);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --eps-up [1,2]",
          "Epsilon scale-up after Metropolis acceptance", opts_.hmc.epsilon_up);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --eps-down [0,1]",
          "Epsilon scale-down after Metropolis rejection", opts_.hmc.epsilon_down);
#ifdef PARALLEL
   fprintf(out_, "  %-30s %s\n", "    --no-sort",
          "Use conventional instead of strictly sorted parallel tempering");
#endif
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
         "Weight of central profile column in CRF initialization", opts_.weight_center);
  fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
          "Parameter for exponential decay of window weights in CRF initialization", 
          opts_.weight_decay);
}

template<class Abc>
void CSHmcApp<Abc>::InitSubstitutionMatrix() {
  BlosumType type = BlosumTypeFromString(opts_.blosum_type);
  sm_.reset(new BlosumMatrix(type));
}

template<class Abc>
void CSHmcApp<Abc>::ReadTrainingData() {
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
void CSHmcApp<Abc>::InitCrf() {
  // Read CRF from jumpstart file?
  if (opts_.modelfile.find(".crf") != std::string::npos) {
    fprintf(out_, "Reading CRF from %s ...\n",
            GetBasename(opts_.modelfile).c_str());

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
    GaussianCrfInit<Abc> init(opts_.gauss_init, *sm_, opts_.hmc.seed);
    crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));

  } else {  // init by sampling from training set
    fputs("Initializing CRF by sampling from training set ...\n", out_);
    MatrixPseudocounts<Abc> pc(*sm_);
    ConstantAdmix admix(opts_.pc_init);
    SamplingCrfInit<Abc, TrainingSequence<Abc> > init(trainset_, pc, admix, *sm_, opts_.hmc.seed, 
        opts_.weight_center, opts_.weight_decay);
    crf_.reset(new Crf<Abc>(opts_.nstates, trainset_.front().x.length(), init));
  }
}

template<class Abc>
int CSHmcApp<Abc>::Run() {
  int myrank = 0;
  InitSubstitutionMatrix();
  ReadTrainingData();
  InitCrf();

#ifdef PARALLEL
  int nreplicas;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nreplicas);

  if (nreplicas > 1) {
    // Append r1, r2, etc. to outfile when running parallel tempering
    if (!opts_.outfile.empty()) {
      string replica = strprintf(".r%d", myrank + 1);
      opts_.outfile.append(replica);
      opts_.crffile.append(replica);
      LOG(ERROR) << opts_.outfile;
      LOG(ERROR) << opts_.crffile;

      fprintf(out_, "Starting replica %d out of %d writing to %s ...\n\n",
              myrank + 1, nreplicas, opts_.outfile.c_str());
    } else {
      throw Exception("Empty outfile when running parallel tempering!");
    }

    // Setup parallel tempering environment if needed
    if (!opts_.no_sort)
      partemp_.reset(new SortedPT(myrank, nreplicas));
    else
      partemp_.reset(new SingleExchPT(myrank, nreplicas, opts_.theta,
                                      opts_.hmc.seed));
  }
#endif

  FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), "w");

  // Setup function objects on training set and validation set
  scoped_ptr<DerivCrfFuncPrior<Abc> > prior;
  if (opts_.prior == 1)
      prior.reset(new GaussianDerivCrfFuncPrior<Abc>(
          opts_.sigma_context, 
          opts_.sigma_decay, 
          opts_.sigma_bias));
  else
      prior.reset(new LassoDerivCrfFuncPrior<Abc>(
          opts_.sigma_context, 
          opts_.sigma_decay, 
          opts_.sigma_bias));
  Likelihood loglike(valset_, *sm_);
  Gradient gradient(trainset_, *sm_, *prior);
  HmcState<Abc> state(*crf_);

  // Setup leapfrog proposal functor
  if (opts_.hmc.sgd_epochs > 0)
    propose_.reset(new BasinHopping(gradient, opts_.hmc, opts_.sgd, myrank));
  else
    propose_.reset(new Leapfrog(gradient, opts_.hmc, myrank));

  // Run 'nsteps' HMC sampling steps
  HmcStep(opts_.nsteps,    // number of HMC steps
          state,           // current sampling state
          *propose_,       // proposal functor
          loglike,         // likelihood functor on validation set
          opts_.hmc.seed,  // seed for Metropolis decision
          partemp_.get(),  // parallel tempering encapsulation
          fout);           // output stream for progress table

  if (!opts_.outfile.empty()) fclose(fout);

  fout = fopen(opts_.crffile.c_str(), "w");
  if (!fout) throw Exception("Can't write to file '%s'!", opts_.crffile.c_str());
  state.crf.Write(fout);
  fclose(fout);
  fprintf(out_, "\nWrote best validation set CRF to %s\n", opts_.crffile.c_str());

  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
#endif
  int rv = cs::CSHmcApp<cs::AA>().main(argc, argv, stdout, "cshmc");
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return rv;
}
