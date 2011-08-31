// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "alignment-inl.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"
#include "crf-inl.h"
#include "getopt_pp.h"
#include "progress_bar.h"
#include "matrix_pseudocounts-inl.h"
#include "sequence-inl.h"
#include "tamura_nei_matrix.h"
#include "training_sequence.h"
#include "training_profile.h"

using namespace GetOpt;
using std::string;
using std::vector;
using std::count;

namespace cs {

struct CSTrainSetAppOptions {
  CSTrainSetAppOptions() { Init(); }

  void Init() {
    dir             = "./";
    dir_col         = "";
    nsamples        = 1000000;
    wlen            = 0;
    seed            = 0;
    pc_admix        = 0.01;
    pc_ali          = 12.0;
    neff_col        = 0.0;
    min_neff        = 0.0;
    min_neff_col    = 0.0;
    max_neff_col    = 20.0;
    max_win         = 1000;
    sampling_mode   = 0;
    weight_center   = 1.6;
    weight_decay    = 0.85;
    pc_engine       = "auto";
    singletons      = 0;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if ((sampling_mode == 1 || sampling_mode == 2 || sampling_mode == 3) 
        && wlen == 0) throw Exception("No window length given!");
    if (outfile.empty()) throw Exception("No output file provided!");
    if (wlen > 0 && !(wlen & 1)) throw Exception("Window length must be odd!");
    if (pc_admix < 0 || pc_admix > 1.0) throw Exception("Pseudocounts admix invalid!");
    if (min_neff < 1.0 && min_neff != 0.0) throw Exception("Minimum Neff invalid!");
    if (min_neff_col < 1.0 && min_neff_col != 0.0) throw Exception("Minimum Neff in pseuocounts column invalid!");
    if (max_neff_col < 1.0 && max_neff_col != 0.0) throw Exception("Maximum Neff in pseuocounts column invalid!");
    if (neff_col < 1.0 && neff_col != 0.0) throw Exception("Neff in pseuocounts column invalid!");
  }

  void PrintOptions(FILE* out) const {
    fprintf(out, "  %-20s: %s\n", "-d, --dir", dir.c_str()); 
    fprintf(out, "  %-20s: %s\n", "-e, --dir-col", dir_col.c_str()); 
    fprintf(out, "  %-20s: %s\n", "-o, --outfile", outfile.c_str()); 
    fprintf(out, "  %-20s: %zu\n", "-N, --size", nsamples); 
    fprintf(out, "  %-20s: %zu\n", "-W, --wlen", wlen); 
    fprintf(out, "  %-20s: %i\n", "-s, --sampling-mode", sampling_mode); 
    fprintf(out, "  %-20s: %s\n", "-D, --context-data", modelfile.c_str()); 
    fprintf(out, "  %-20s: %.2f\n", "-x, --pc-admix", pc_admix); 
    fprintf(out, "  %-20s: %.1f\n", "-c, --pc-ali", pc_ali); 
    fprintf(out, "  %-20s: %.1f\n", "-y, --neff-col", neff_col); 
    fprintf(out, "  %-20s: %.1f\n", "-n, --min-neff", min_neff); 
    fprintf(out, "  %-20s: %.1f\n", "-m, --min-neff-col",min_neff_col); 
    fprintf(out, "  %-20s: %.1f\n", "-M, --max-neff-col",max_neff_col); 
    fprintf(out, "  %-20s: %i\n", "-r, --seed", seed); 
    fprintf(out, "  %-20s: %.2f\n", "-g, --singletons", singletons); 
  }
 
  string dir;            // directory with input data to build trainset from
  string dir_col;         // directory with count profiles for training pseudocounts column
  string ext;            // file extension of input files with counts
  string outfile;        // created training set
  string modelfile;      // input file with context profile library or HMM
  string pc_engine;      // pseudocount engine
  size_t nsamples;       // size of training set to create
  size_t wlen;           // length of context window (zero means full-length)
  double pc_admix;       // minimal pseudocount admix for count profiles
  double pc_ali;         // constant in pseudocount calculation for alignments
  double neff_col;       // target Neff in training pseudocounts column
  double min_neff;       // minimum Neff in count profiles
  double min_neff_col;   // minimum Neff in in training pseudocounts column
  double max_neff_col;   // maximum Neff in in training pseudocounts column
  size_t max_win;        // maximal number of windows per protein
  int sampling_mode;     // sample sequence-column pairs from alignments
  double weight_center;  // weight of central column in multinomial emission
  double weight_decay;   // exponential decay of window weights
  unsigned int seed;     // seed
  double singletons;     // fraction of singletons in training profiles
};  // CSTrainSetAppOptions


template<class Abc>
class CSTrainSetApp : public Application {
 protected:
  typedef vector<CountProfile<Abc> > ProfileSet;
  typedef vector<TrainingSequence<Abc> > TrainSeqs;
  typedef vector<TrainingProfile<Abc> > TrainProfiles;
  typedef vector<string> Files;

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
  // Samples full-length profiles from database of profiles.
  void SampleFullLengthProfiles(ProfileSet& samples);
  // Samples profile windows from database of profiles.
  void SampleProfileWindows(ProfileSet& samples);
  // Samples sequence-column pairs from database of profiles and alignments.
  void SampleTrainingSeqs(TrainSeqs& samples, size_t nsamples_);
  // Samples profile-column pairs from database of profiles and alignments.
  void SampleTrainingProfiles(TrainProfiles& samples);
  // Initializes substitution matrix (specialized by alphabet type).
  void InitSubstitutionMatrix();

  CSTrainSetAppOptions opts_;
  ProfileSet profiles_;
  ProfileSet* profiles_col_;
  Files files_;
  scoped_ptr<SubstitutionMatrix<Abc> > sm_;
  scoped_ptr<ContextLibrary<Abc> > lib_;
  scoped_ptr<Crf<Abc> > crf_;
  scoped_ptr<Pseudocounts<Abc> > pc_;
  vector<vector<size_t> > stats_;

};  // CSTrainSetApp


template<class Abc>
void CSTrainSetApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('d', "dir", opts_.dir, opts_.dir);
  ops >> Option('e', "dir-col", opts_.dir_col, opts_.dir_col);
  ops >> Option('N', "size", opts_.nsamples, opts_.nsamples);
  ops >> Option('W', "wlen", opts_.wlen, opts_.wlen);
  ops >> Option('r', "seed", opts_.seed, opts_.seed);
  ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  ops >> Option('y', "neff-col", opts_.neff_col, opts_.neff_col);
  ops >> Option('n', "min-neff", opts_.min_neff, opts_.min_neff);
  ops >> Option('m', "min-neff-col", opts_.min_neff_col, opts_.min_neff_col);
  ops >> Option('M', "max-neff-col", opts_.max_neff_col, opts_.max_neff_col);
  ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
  ops >> Option('s', "sampling-mode", opts_.sampling_mode, opts_.sampling_mode);
  ops >> Option('g', "singletons", opts_.singletons, opts_.singletons);

  opts_.Validate();

  if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
    opts_.pc_engine = GetFileExt(opts_.modelfile);
}

template<class Abc>
void CSTrainSetApp<Abc>::PrintBanner() const {
  fputs("Build a training set by sampling from a larger database.\n",
        out_);
}

template<class Abc>
void CSTrainSetApp<Abc>::PrintUsage() const {
  fputs("Usage: cstrainset -o <outfile> [options]\n", out_);
}

template<class Abc>
void CSTrainSetApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with sampled training set");
  fprintf(out_, "  %-30s %s (def=%s)\n", "-d, --dir <dir>",
          "Directory with count profiles (and A3M alignments)", opts_.dir.c_str());
  fprintf(out_, "  %-30s %s (def=%s)\n", "-e, --dir-col <dir>",
          "Directory with count profiles for training pseudocounts column", opts_.dir_col.c_str());
  fprintf(out_, "  %-30s %s (def=%zu)\n", "-N, --size [0,inf[",
          "Size of training set to sample", opts_.nsamples);
  fprintf(out_, "  %-30s %s\n", "-W, --wlen [0,inf[",
          "Sample profile windows of length W instead of full-length");
  fprintf(out_, "  %-30s %s (def=%d)\n", "-s, --sampling-mode [0-3]",
          "Sampling mode", opts_.sampling_mode);
  fprintf(out_, "  %-30s %s\n", "", "0: sample count profiles");
  fprintf(out_, "  %-30s %s\n", "", "1: sample training seqs from alignment query");
  fprintf(out_, "  %-30s %s\n", "", "2: sample training seqs from full alignment");
  fprintf(out_, "  %-30s %s\n", "", "3: sample training profiles");
  fprintf(out_, "  %-30s %s\n", "-D, --context-data <file>",
          "Context data for adding cs-pseudocounts instead of subst. matrix PCs");
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix to be added to count profiles", opts_.pc_admix);
  fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant in pseudocount calculation for alignments", opts_.pc_ali);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-y, --neff-col [1,inf[",
          "Target number of effective sequences in training pseudocounts column", opts_.neff_col);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-n, --min-neff [1,20]",
          "Minimum number of effective sequences in profiles", opts_.min_neff);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-m, --min-neff-col [1,20]",
          "Minimum number of effective sequences in training pseudocounts column", opts_.min_neff_col);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-M, --max-neff-col [1,20]",
          "Maximum number of effective sequences in training pseudocounts column", opts_.max_neff_col);
  fprintf(out_, "  %-30s %s (def=%u)\n", "-r, --seed [0,inf[",
          "Seed for random number generator", opts_.seed);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-g, --singletons [0,1]",
          "Fraction of singletons in training profiles", opts_.singletons);
}

template<class Abc>
void CSTrainSetApp<Abc>::InitSubstitutionMatrix() {
  sm_.reset(new BlosumMatrix());
}

template<>
void CSTrainSetApp<Dna>::InitSubstitutionMatrix() {
  sm_.reset(new TamuraNeiMatrix());
}

template<class Abc>
int CSTrainSetApp<Abc>::Run() {
  fputs("Running cstrainset with following parameters:\n", out_);
  opts_.PrintOptions(out_);
  fputs("\n", out_);

  InitSubstitutionMatrix();
  Ran ran(opts_.seed);  // random number generator for sorting

  // Setup pseudocount engine
  if (!opts_.modelfile.empty() && opts_.pc_engine == "lib") {
    fprintf(out_, "Reading context library from %s ...\n",
            GetBasename(opts_.modelfile).c_str());
    FILE* fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    lib_.reset(new ContextLibrary<Abc>(fin));
    TransformToLog(*lib_);
    fclose(fin);
    pc_.reset(new LibraryPseudocounts<Abc>(*lib_, opts_.weight_center,
                                              opts_.weight_decay));

  } else if (!opts_.modelfile.empty() && opts_.pc_engine == "crf") {
    fprintf(out_, "Reading CRF from %s ...\n",
            GetBasename(opts_.modelfile).c_str());
    FILE* fin = fopen(opts_.modelfile.c_str(), "r");
    if (!fin)
      throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
    crf_.reset(new Crf<Abc>(fin));
    fclose(fin);
    pc_.reset(new CrfPseudocounts<Abc>(*crf_));
  } else {
    pc_.reset(new MatrixPseudocounts<Abc>(*sm_));
  }

  // Glob count-profile files with .prf extension wich are both in 
  // opts_.dir and opts_.dir_col
  Files files_glob;
  fputs("Globbing input directory for profiles ...", out_);
  fflush(out_);
  GetAllFiles(opts_.dir, files_glob, "prf");  
  random_shuffle(files_glob.begin(), files_glob.end(), ran);
  if (!opts_.dir_col.empty()) {
    Files files;
    std::map<string, bool> files_col;
    GetAllFiles(opts_.dir_col, files, "prf");
    for (Files::iterator it = files.begin(); it != files.end(); ++it) 
      files_col[*it] = true;
    files.clear();
    for (Files::iterator it = files_glob.begin(); it != files_glob.end(); ++it) {
      if (files_col.count(*it) > 0) files.push_back(*it);
    }
    files_glob = files;
  }
  fprintf(out_, "\n%zu files globbed\n", files_glob.size());

  // Read all count profiles which fulfil filter criteria
  fprintf(out_, "Reading %zu profiles into memory ...\n", files_glob.size());
  ProgressBar progress(out_, 70, files_glob.size());
  if (opts_.dir_col.empty()) 
    profiles_col_ = &profiles_; // Use a single set of count profiles
  else
    profiles_col_ = new ProfileSet(); // Use two separate sets of counts profiles
  for (Files::iterator it = files_glob.begin(); it != files_glob.end(); ++it) {
    // Read next count profile from opts_.dir
    string filename = opts_.dir + kDirSep + *it;
    FILE* fin = fopen(filename.c_str(), "r");
    if (!fin) throw Exception("Unable to open file '%s'!", filename.c_str());
    CountProfile<Abc> prof(fin);
    fclose(fin);
    progress.Advance();

    // Discard count profiles which are too short or with Neff < opts_.min_neff
    double neff = Neff(prof);
    if (neff < opts_.min_neff || prof.counts.length() < opts_.wlen)
      continue;

    if (opts_.dir_col.empty()) {
      // Discard count profiles with Neff < opts_.min_neff_col
      if (neff < opts_.min_neff_col || neff > opts_.max_neff_col) continue;
      files_.push_back(*it);
      profiles_.push_back(prof);

    } else {
      // Read next count profile from opts_.dir_col
      filename = opts_.dir_col + kDirSep + *it;
      fin = fopen(filename.c_str(), "r");
      if (!fin) throw Exception("Unable to open file '%s'!", filename.c_str());
      CountProfile<Abc> prof_col(fin);
      fclose(fin);
      if (prof_col.counts.length() != prof.counts.length()) {
        throw Exception("Size of profile '%s' does not match the size of the corresponding profile!",
            filename.c_str());
      }
      // Discard count profiles with Neff < opts_.min_neff_col
      double neff_col = Neff(prof_col);
      if (neff_col < opts_.min_neff_col || neff_col > opts_.max_neff_col) continue;
      files_.push_back(*it);
      profiles_.push_back(prof);
      profiles_col_->push_back(prof_col);
    }
  }
  for (size_t i = 0; i < profiles_.size(); ++i)
    stats_.push_back(vector<size_t>(profiles_[i].length(), 0));
  assert(files_.size() == profiles_.size());
  assert(profiles_.size() == profiles_col_->size());
  fprintf(out_, "\n%zu profiles passed filter\n", profiles_.size());

  // Add pseudocounts to count profiles in parallel
  if (opts_.pc_admix > 0.0) {
    fputs("Adding pseudocounts to profiles ...\n", out_);
    CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
    const int nprof = profiles_.size();
    progress.Init(nprof);
#pragma omp parallel for schedule(static)
    for (int n = 0; n < nprof; ++n) {
      profiles_[n].counts = pc_->AddTo(profiles_[n], admix);
      Normalize(profiles_[n].counts, profiles_[n].neff);
      if (!opts_.dir_col.empty()) {
        CountProfile<Abc>& prof = profiles_col_->at(n);
        prof.counts = pc_->AddTo(prof, admix);
        Normalize(prof.counts, prof.neff);
      }
#pragma omp critical (pseudocounts_advance_progress)
      progress.Advance();
    }
    fputs("\n", out_);
  }

  // Now we can call the actual sampling method
  if (opts_.sampling_mode == 0) {
    ProfileSet samples;
    samples.reserve(opts_.nsamples);
    if (opts_.wlen == 0)
      SampleFullLengthProfiles(samples);
    else
      SampleProfileWindows(samples);
    FILE* fout = fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
    fprintf(out_, "Writing %zu count profiles to %s ...", samples.size(), opts_.outfile.c_str());
    fflush(out_);

  } else if (opts_.sampling_mode == 1 || opts_.sampling_mode == 2) {
    TrainSeqs samples;
    samples.reserve(opts_.nsamples);
    SampleTrainingSeqs(samples, opts_.nsamples);
    FILE* fout = fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
    fprintf(out_, "Writing %zu training sequences to %s ...", samples.size(),
        opts_.outfile.c_str());
    fflush(out_);
    WriteAll(samples, fout);
    fclose(fout);

  } else {
    TrainProfiles samples;
    samples.reserve(opts_.nsamples);
    SampleTrainingProfiles(samples);
    FILE* fout = fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
    fprintf(out_, "Writing %zu training profiles to %s ...", samples.size(),
        opts_.outfile.c_str());
    fflush(out_);
    WriteAll(samples, fout);
    fclose(fout);

  }
  // Calculate sampling statistics
  if (stats_.size() > 0) {
    double stats_profiles[] = {0.0, 0.0, 0.0, DBL_MAX, -DBL_MAX};
    double stats_cols[] = {0.0, 0.0, 0.0, DBL_MAX, -DBL_MAX};
    size_t nused_cols = 0;
    for (size_t n = 0; n < stats_.size(); ++n) {
      size_t nused = 0;
      size_t nsamples = 0;
      for (size_t m = 0; m < stats_[n].size(); ++m) {
        if (stats_[n][m] > 0) {
          nused++;
          nsamples += stats_[n][m];
          stats_cols[1] += stats_[n][m];
          stats_cols[3] = MIN(stats_cols[3], stats_[n][m]);
          stats_cols[4] = MAX(stats_cols[4], stats_[n][m]);
        }
      }
      stats_cols[0] += static_cast<double>(nused) / stats_[n].size();
      nused_cols += nused;
      if (nsamples > 0) {
        stats_profiles[0]++;
        stats_profiles[1] += nsamples;
        stats_profiles[3] = MIN(stats_profiles[3], nsamples);
        stats_profiles[4] = MAX(stats_profiles[4], nsamples);
      }
    }
    stats_cols[0] /= stats_.size();
    stats_cols[1] /= nused_cols;
    stats_profiles[1] /= stats_profiles[0]; 
    for (size_t n = 0; n < stats_.size(); ++n) {
      size_t nsamples = 0;
      for (size_t m = 0; m < stats_[n].size(); ++m) {
        if (stats_[n][m] > 0) {
          nsamples += stats_[n][m];
          stats_cols[2] += pow(stats_[n][m] - stats_cols[1], 2);
        }
      }
      if (nsamples > 0)
        stats_profiles[2] += pow(nsamples - stats_profiles[1], 2);
    }
    stats_cols[2] = sqrt(stats_cols[2] / nused_cols);
    stats_profiles[2] = sqrt(stats_profiles[2] / stats_profiles[0]);
    stats_profiles[0] /= stats_.size();
    fprintf(out_, "\nSampling statistics:\n");
    fprintf(out_, " %-20s %8s %8s %8s %8s %8s\n", "", "used", "avg", "std", "min", "max");
    fprintf(out_, " %-20s %8.2f %8.2f %8.2f %8d %8d\n", "Samples per profile", 
        stats_profiles[0], stats_profiles[1], stats_profiles[2], 
        static_cast<int>(stats_profiles[3]), static_cast<int>(stats_profiles[4]));
    fprintf(out_, " %-20s %8.2f %8.2f %8.2f %8d %8d\n", "Samples per column", 
        stats_cols[0], stats_cols[1], stats_cols[2], 
        static_cast<int>(stats_cols[3]), static_cast<int>(stats_cols[4]));
  }

  if (opts_.dir_col.length()) delete profiles_col_;
  fputs("\nDone!\n", out_);

  return 0;
}


template<class Abc>
void CSTrainSetApp<Abc>::SampleFullLengthProfiles(ProfileSet& samples) {
  const size_t nsamples = MIN(profiles_.size(), opts_.nsamples);
  fprintf(out_, "Sampling %zu profiles ...\n", nsamples);
  ProgressBar progress(out_, 70, nsamples);

  for (size_t n = 0; n < profiles_.size(); ++n) {
    if (samples.size() < nsamples) {
      samples.push_back(profiles_[n]);
      progress.Advance();
    }
  }

  fputs("\nShuffling sampled profiles set ...", out_);
  fflush(out_);
  Ran ran(opts_.seed);
  random_shuffle(samples.begin(), samples.end(), ran);
  fputs("\n", out_);
}

template<class Abc>
void CSTrainSetApp<Abc>::SampleProfileWindows(ProfileSet& samples) {
  // Precompute total pool size taking masking into account
  const size_t center = (opts_.wlen - 1) / 2;
  size_t nall = 0;
  for (size_t n = 0; n < profiles_.size(); ++n) {
    const CountProfile<Abc>& cp = profiles_[n];
    const size_t ncol =  cp.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (cp.neff[i + center] >= opts_.min_neff) unmasked++;
    }
    nall += MIN(opts_.max_win, unmasked);
  }
  const size_t nsamples = MIN(opts_.nsamples, nall);

  fprintf(out_, "Sampling %zu profiles with W=%zu out of %zu windows ...\n",
          nsamples, opts_.wlen, nall);

  // Iterate over input data and build counts profiles
  ProgressBar progress(out_, 70, nsamples);
#pragma omp parallel for schedule(static)
  for (size_t n = 0; n < profiles_.size(); ++n) {
    if (samples.size() < nsamples) {
      Ran ran(opts_.seed + n);
      const CountProfile<Abc>& cp = profiles_[n];
      LOG(ERROR) << "sampling " << files_[n];

      vector<int> shuffle;
      for (size_t i = 0; i <= cp.counts.length() - opts_.wlen; ++i)
        if (cp.neff[i + center] >= opts_.min_neff) shuffle.push_back(i);
      random_shuffle(shuffle.begin(), shuffle.end(), ran);

      // Proportional number of training windows we need to sample
      size_t w = MIN(opts_.max_win, shuffle.size());
      size_t s = ceil(nsamples * (static_cast<double>(w) / nall));
      LOG(ERROR) << strprintf("n=%zu nsamples=%zu  nall=%zu  len=%zu  s=%zu",
                              n, nsamples, nall, cp.counts.length(), s);

      // Sample only a fraction of the profile indices.
      shuffle.erase(shuffle.begin() + s, shuffle.end());
      // Copy profile windows into 'samples' vector
      for (size_t i = 0; i < shuffle.size(); ++i) {
        CountProfile<Abc> p(cp, shuffle[i], opts_.wlen);
#pragma omp critical (add_sample)
        if (samples.size() < nsamples) {
          samples.push_back(p);
          stats_[n][shuffle[i]]++;
          progress.Advance();
        }
      }
    }
  }

  fputs("\nShuffling sampled profiles set ...", out_);
  fflush(out_);
  Ran ran(opts_.seed);
  random_shuffle(samples.begin(), samples.end(), ran);
  fputs("\n", out_);
}

template<class Abc>
void CSTrainSetApp<Abc>::SampleTrainingSeqs(TrainSeqs& samples, size_t nsamples_) {
  const size_t center = (opts_.wlen - 1) / 2;
  Vector<double> neff(profiles_.size(), 0.0);
  double nall = 0.0;  // total pool size

  // Precompute Neff and total pool size taking masking into account
  for (size_t n = 0; n < profiles_.size(); ++n) {
    const CountProfile<Abc>& cp = profiles_[n];
    const size_t ncol = cp.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (cp.neff[i + center] >= opts_.min_neff &&
          profiles_col_->at(n).neff[i + center] >= opts_.min_neff_col &&
          profiles_col_->at(n).neff[i + center] <= opts_.max_neff_col) {
        unmasked++;
        neff[n] += opts_.sampling_mode == 2 ? cp.neff[i + center] : 1.0;
      }
    }
    if (unmasked > 0) 
      neff[n] /= unmasked;  // average Neff in columns that passed Neff-filter
    nall += MIN(opts_.max_win, unmasked) * neff[n];
  }
  const size_t nsamples = MIN(nsamples_, floor(nall));

  fprintf(out_, "Sampling %zu training sequences out of %.0f windows ...\n",
          nsamples, nall);

  // Iterate over profile-alignment pairs and sample from each pair a certain
  // number of sequence windows together with the corresponding counts column
  // at the central window position.
  ProgressBar progress(out_, 70, nsamples);
#pragma omp parallel for schedule(static)
  for (size_t n = 0; n < profiles_.size(); ++n) {
    if (samples.size() < nsamples) {
      Ran ran(opts_.seed + n);
      const CountProfile<Abc>& cp = profiles_[n];
      size_t nmatch = cp.counts.length();

      // Filter out columns that don't suffice Neff-filter, same as above
      size_t ncol = nmatch - opts_.wlen + 1;
      Vector<bool> masked(ncol, false);
      int unmasked = 0;
      for (size_t i = 0; i < ncol; ++i) {
        if (cp.neff[i + center] < opts_.min_neff ||
            profiles_col_->at(n).neff[i + center] < opts_.min_neff_col ||
            profiles_col_->at(n).neff[i + center] > opts_.max_neff_col) masked[i] = true;
        else unmasked++;
      }

      // Determine workload for this profile
      double nprf = MIN(opts_.max_win, static_cast<size_t>(unmasked)) * neff[n];
      double todo = MIN(1.0, nprf / MAX(0.0, nall)) * (nsamples - samples.size());
      double left = unmasked * neff[n];

      if (todo == 0.0 || samples.size() >= nsamples) continue;

      LOG(ERROR) <<
        strprintf("n=%4zu neff=%4.1f L=%4zu nall=%9.1f nprf=%9.1f todo=%5.1f",
                  n, neff[n], nmatch, nall, nprf, todo);

#pragma omp atomic
      nall -= nprf;

      // Read-in the underlying A3M alignment of profile 'n'
      string file = opts_.dir + kDirSep + GetBasename(files_[n], false) + ".a3m";
      FILE* fp = fopen(file.c_str(), "r");
      if (!fp) throw Exception("Can't open alignment file '%s'!", file.c_str());
      Alignment<Abc> ali(fp, A3M_ALIGNMENT);
      fclose(fp);

      // Check if number of match columns in profile and alignment are correct
      if (nmatch != ali.nmatch())
        throw Exception("Number of matchcols in ali '%s' should be %zu but is %zu!",
                        GetBasename(file).c_str(), cp.counts.length(), nmatch);

      // Try out sequences within each column in random order
      vector<int> shuffle;
      if (opts_.sampling_mode == 2) {
        for (size_t k = 0; k < ali.nseqs(); ++k) 
          shuffle.push_back(k);
      } else {
        shuffle.push_back(0);
      }

      while (todo > 0.0 && samples.size() < nsamples && unmasked > 0) {
        // Pick a random column as window start position
        size_t m = ran(ncol);
        if (masked[m]) continue;
        // Determine sampling work load for column 'm' (10% more as safe-guard)
        double neff_center = opts_.sampling_mode == 2 ? cp.neff[m + center] : 1.0;
        double s = 1.1 * todo * (neff_center / left);
        LOG(ERROR) << strprintf("todo=%5.2f  neff[%2zu]=%4.1f  left=%9.1f  s=%4.2f",
                                todo, m, cp.neff[m + center], left, s);
        assert(s > 0.0);
        // Shuffle sequences in column 'm'
        random_shuffle(shuffle.begin(), shuffle.end(), ran);
        size_t l = 0;  // index in shuffle vector

        while (ran.doub() < s && l < shuffle.size() && samples.size() < nsamples) {
          if (!masked[m]) {
            masked[m] = true;
            unmasked--;
            left -= neff_center;
          }

          while (l < shuffle.size()) {
            size_t k = shuffle[l++];
            LOG(ERROR) << strprintf("m=%zu: trying k=%zu", m, k);

            // Build sequence window along with checking if it's valid
            Sequence<Abc> seq(opts_.wlen);
            string hdr(ali.header(k));
            string::size_type si = hdr.find(' ');
            seq.set_header(hdr.substr(0, si != string::npos ? si : hdr.size()));
            bool valid = true;
            size_t i = ali.col_idx(m);  // global column index of match column 'm'
            size_t j = 0;
            while (j < opts_.wlen && i < ali.ncols()) {
              if ((ali.is_match(i) && ali(k,i) >= Abc::kGap) ||
                  (!ali.is_match(i) && ali(k,i) < Abc::kGap)) {
                valid = false;
                break;
              }
              if (ali.is_match(i)) seq[j++] = ali(k,i);
              i++;
            }

            if (valid && samples.size() < nsamples) {
              assert(j == opts_.wlen);
              // Cut training pseudocounts column from profiles_col_
              CountProfile<Abc> cp_col = CountProfile<Abc>(profiles_col_->at(n), m, opts_.wlen);
              if (opts_.neff_col != 0.0 && cp_col.neff[center] < opts_.neff_col) {               
                Profile<Abc> p = cp_col.counts;
                Normalize(p, 1.0);
                // Estimated target Neff in count profile = Neff / Neff(cp_i) * Neff(p)
                cp_col.counts = pc_->AddTo(cp_col, opts_.neff_col / cp_col.neff[center] * Neff(p));
                Assign(cp_col.neff, Neff(cp_col.counts));
                Normalize(cp_col.counts, cp_col.neff);
              }
              ProfileColumn<Abc> col(cp_col.counts[center]);
              if (samples.size() < nsamples) {
#pragma omp critical (add_sample)
                {
                  samples.push_back(TrainingSequence<Abc>(seq, col));
                  stats_[n][m]++;
                  progress.Advance();
                }
                s -= 1;     // reduce work load for this column by one
                todo -= 1;  // reduce work load for whole profile by one
                LOG(ERROR) << strprintf("SUCCESS!  todo=%5.2f  s=%4.2f seqs=%zu",
                                        todo, s, shuffle.size() - l);
                break;      // pick next sequence window
              }
            }
          }
        }
      }
    }
  }

  fputs("\nShuffling sampled training sequences ...", out_);
  fflush(out_);
  Ran ran(opts_.seed);
  random_shuffle(samples.begin(), samples.end(), ran);
  fputs("\n", out_);
}

template<class Abc>
void CSTrainSetApp<Abc>::SampleTrainingProfiles(TrainProfiles& samples) {

  size_t nsingletons = static_cast<size_t>(opts_.singletons * opts_.nsamples);
  if (nsingletons > 0) {
    // Sample singletons
    TrainSeqs samples_seqs;
    samples_seqs.reserve(nsingletons);
    SampleTrainingSeqs(samples_seqs, nsingletons);
    nsingletons = samples_seqs.size();
    for (size_t n = 0; n < nsingletons; ++n)
      samples.push_back(TrainingProfile<Abc>(samples_seqs[n].x, samples_seqs[n].y));
  }

  // Precompute total pool size taking masking into account
  const size_t center = (opts_.wlen - 1) / 2;
  size_t nall = 0;
  for (size_t n = 0; n < profiles_.size(); ++n) {
    const CountProfile<Abc>& cp = profiles_[n];
    const size_t ncol =  cp.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (cp.neff[i + center] >= opts_.min_neff &&
          profiles_col_->at(n).neff[i + center] >= opts_.min_neff_col &&
          profiles_col_->at(n).neff[i + center] <= opts_.max_neff_col) unmasked++;
    }
    nall += MIN(opts_.max_win, unmasked);
  }
  const size_t nsamples = MIN(opts_.nsamples - nsingletons, nall);

  fprintf(out_, "Sampling %zu training profiles with W=%zu out of %zu windows ...\n",
          nsamples, opts_.wlen, nall);

  // Iterate over input data and build training profiles
  ProgressBar progress(out_, 70, nsamples);
#pragma omp parallel for schedule(static)
  for (size_t n = 0; n < profiles_.size(); ++n) {
    if (samples.size() < nsamples + nsingletons) {
      Ran ran(opts_.seed + n + 1);
      const CountProfile<Abc>& cp = profiles_[n];
      LOG(ERROR) << "sampling " << files_[n];

      vector<int> shuffle;
      for (size_t i = 0; i <= cp.counts.length() - opts_.wlen; ++i) {
        if (cp.neff[i + center] >= opts_.min_neff &&
            profiles_col_->at(n).neff[i + center] >= opts_.min_neff_col &&
            profiles_col_->at(n).neff[i + center] <= opts_.max_neff_col) shuffle.push_back(i);
      }
      random_shuffle(shuffle.begin(), shuffle.end(), ran);

      // Proportional number of training windows we need to sample
      size_t w = MIN(opts_.max_win, shuffle.size());
      size_t s = ceil(nsamples * (static_cast<double>(w) / nall));
      LOG(ERROR) << strprintf("n=%zu nsamples=%zu  nall=%zu  len=%zu  s=%zu",
                              n, nsamples, nall, cp.counts.length(), s);
      s = MIN(s, shuffle.size());

      // Copy profile windows into 'samples' vector
      for (size_t i = 0; i < s; ++i) {
        CountProfile<Abc> p(cp, shuffle[i], opts_.wlen);
        // Cut training pseudocounts column from profiles_col_
        CountProfile<Abc> cp_col = CountProfile<Abc>(profiles_col_->at(n), shuffle[i], opts_.wlen);
        if (opts_.neff_col != 0.0 && cp_col.neff[center] < opts_.neff_col) {               
          Profile<Abc> p = cp_col.counts;
          Normalize(p, 1.0);
          // Estimated target Neff in count profile = Neff / Neff(M_i) * Neff(CP_i)
          p = pc_->AddTo(cp_col, opts_.neff_col / cp_col.neff[center] * Neff(p));
          cp_col.counts = p;
          Assign(cp_col.neff, Neff(p));
          Normalize(cp_col.counts, cp_col.neff);
        }
        ProfileColumn<Abc> col(cp_col.counts[center]);
#pragma omp critical (add_sample)
        if (samples.size() < nsamples + nsingletons) {
          Normalize(p.counts, p.neff);
          samples.push_back(TrainingProfile<Abc>(p, col));
          stats_[n][shuffle[i]]++;
          progress.Advance();
        }
      }
    }
  }

  fputs("\nShuffling sampled training profiles ...", out_);
  fflush(out_);
  Ran ran(opts_.seed);
  random_shuffle(samples.begin(), samples.end(), ran);
  fputs("\n", out_);
}


}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSTrainSetApp<cs::Dna>().main(argc, argv, stdout, "cstrainset");
  else
    return cs::CSTrainSetApp<cs::AA>().main(argc, argv, stdout, "cstrainset");
}

