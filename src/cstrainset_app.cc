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
    dir           = "./";
    profile_ext   = "prf";
    ali_ext       = "a3m";
    nsamples      = 1000000;
    wlen          = 0;
    seed          = 0;
    pc_admix      = 0.01;
    pc_ali        = 12.0;
    neff_x_min    = kNeffMin;
    neff_x_max    = kNeffMax;
    neff_y_min    = kNeffMin;
    neff_y_max    = kNeffMax;
    neff_d_min    = 0.0;
    neff_d_max    = 0.0;
    neff_y_target = 0.0;
    round         = 0;
    use_min_round = 0;
    max_win       = 1000;
    sampling_mode = 0;
    weight_center = 1.6;
    weight_decay  = 0.85;
    pc_engine     = "auto";
    singletons    = 0;
    central_mod   = false;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if ((sampling_mode == 1 || sampling_mode == 2 || sampling_mode == 3) 
        && wlen == 0) throw Exception("No window length given!");
    if (outfile.empty()) throw Exception("No output file provided!");
    if (wlen > 0 && !(wlen & 1)) throw Exception("Window length must be odd!");
    if (pc_admix < 0 || pc_admix > 1.0) throw Exception("Pseudocounts admix invalid!");
    if (neff_x_min < kNeffMin) throw Exception("Minimum Neff in profiles invalid!");
    if (neff_x_max < kNeffMin) throw Exception("Maximum Neff in profiles invalid!");
    if (neff_y_min < kNeffMin) throw Exception("Minimum Neff in pseuocounts column invalid!");
    if (neff_y_max < kNeffMin) throw Exception("Maximum Neff in pseuocounts column invalid!");
    if (neff_y_target < kNeffMin && neff_y_target != 0.0) 
      throw Exception("Target Neff in pseuocounts column invalid!");
    if (neff_d_min > 0.0 && neff_d_max == 0.0) neff_d_max = kNeffMax;
  }

  void PrintOptions(FILE* out) const {
    fprintf(out, "  %-20s: %s\n", "-d, --dir", dir.c_str()); 
    fprintf(out, "  %-20s: %s\n", "-o, --outfile", outfile.c_str()); 
    fprintf(out, "  %-20s: %zu\n", "-N, --size", nsamples); 
    fprintf(out, "  %-20s: %zu\n", "-W, --wlen", wlen); 
    fprintf(out, "  %-20s: %i\n", "-s, --sampling-mode", sampling_mode); 
    fprintf(out, "  %-20s: %.2f\n", "-g, --singletons", singletons); 
    fprintf(out, "  %-20s: %s\n", "-D, --context-data", modelfile.c_str()); 
    fprintf(out, "  %-20s: %.2f\n", "-x, --pc-admix", pc_admix); 
    fprintf(out, "  %-20s: %.1f\n", "-c, --pc-ali", pc_ali); 
    fprintf(out, "  %-20s: %.1f\n", "-u, --neff-x-min", neff_x_min); 
    fprintf(out, "  %-20s: %.1f\n", "-U, --neff-x-max", neff_x_max); 
    fprintf(out, "  %-20s: %.1f\n", "-v, --neff-y-min", neff_y_min); 
    fprintf(out, "  %-20s: %.1f\n", "-V, --neff-y-max", neff_y_max); 
    fprintf(out, "  %-20s: %.1f\n", "-y, --neff-y-target", neff_y_target);
    fprintf(out, "  %-20s: %.1f\n", "-j, --neff-d-min", neff_d_min); 
    fprintf(out, "  %-20s: %.1f\n", "-J, --neff-d-max", neff_d_max); 
    fprintf(out, "  %-20s: %zu\n", "-R, --round", round);
    fprintf(out, "  %-20s: %d\n", "-M, --use-min-round", use_min_round);
    fprintf(out, "  %-20s: %d\n", "-C, --central-mod", central_mod);
    fprintf(out, "  %-20s: %d\n", "-r, --seed", seed); 
  }
 
  string dir;           // directory with input data to build trainset from
  string profile_ext;   // file extension of profiles
  string ali_ext;       // file extension of alignments
  string outfile;       // created training set
  string modelfile;     // input file with context profile library or HMM
  string pc_engine;     // pseudocount engine
  size_t nsamples;      // size of training set to create
  size_t wlen;          // length of context window (zero means full-length)
  double pc_admix;      // minimal pseudocount admix for count profiles
  double pc_ali;        // constant in pseudocount calculation for alignments
  double neff_x_min;    // minimum Neff in profiles
  double neff_x_max;    // maximum Neff in profiles
  double neff_y_min;    // minimum Neff in training pseudocounts column
  double neff_y_max;    // maximum Neff in training pseudocounts column
  double neff_y_target; // target Neff in training pseudocounts column
  double neff_d_min;    // minimum Neff in training pseudocounts column
  double neff_d_max;    // maximum Neff in training pseudocounts column
  size_t round;         // PSI-BLAST round for filtering input profiles
  bool use_min_round;   // If several profiles meet the filtering condition, use that of the minimum round
  size_t max_win;       // maximal number of windows per protein
  int sampling_mode;    // sample sequence-column pairs from alignments
  double weight_center; // weight of central column in multinomial emission
  double weight_decay;  // exponential decay of window weights
  bool central_mod;     // modify central profile columns
  unsigned int seed;    // seed
  double singletons;    // fraction of singletons in training profiles

  static const double kNeffMin = 1.0;
  static const double kNeffMax = 20.0;
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
  // Reads profiles which meet the filtering condition into memory.
  void ReadProfiles(ProgressBar& progres);
  // Computes sampling statistics.
  void ComputeStatistics(ProgressBar& progress);
  // Samples full-length profiles from database of profiles.
  void SampleFullLengthProfiles(ProfileSet& samples);
  // Samples profile windows from database of profiles.
  void SampleProfileWindows(ProfileSet& samples);
  // Checks whether the column diversities meet the filtering conditions.
  inline bool IsValidPos(double neff_x, double neff_y) const;
  // Samples sequence-column pairs from database of profiles and alignments.
  void SampleTrainingSeqs(TrainSeqs& samples, size_t nsamples_);
  // Samples profile-column pairs from database of profiles and alignments.
  void SampleTrainingProfiles(TrainProfiles& samples);
  // Initializes substitution matrix (specialized by alphabet type).
  void InitSubstitutionMatrix();

  CSTrainSetAppOptions opts_;
  ProfileSet profiles_x_;
  ProfileSet* profiles_y_;
  bool profiles_xy_;    
  Files files_x_;
  scoped_ptr<SubstitutionMatrix<Abc> > sm_;
  scoped_ptr<ContextLibrary<Abc> > lib_;
  scoped_ptr<Crf<Abc> > crf_;
  scoped_ptr<Pseudocounts<Abc> > pc_;
  vector<vector<size_t> > stats_;

};  // CSTrainSetApp


template<class Abc>
void CSTrainSetApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('d', "dir", opts_.dir, opts_.dir);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('N', "size", opts_.nsamples, opts_.nsamples);
  ops >> Option('W', "wlen", opts_.wlen, opts_.wlen);
  ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
  ops >> Option('s', "sampling-mode", opts_.sampling_mode, opts_.sampling_mode);
  ops >> Option('g', "singletons", opts_.singletons, opts_.singletons);
  ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
  ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
  ops >> Option('u', "neff-x-min", opts_.neff_x_min, opts_.neff_x_min);
  ops >> Option('U', "neff-x-max", opts_.neff_x_max, opts_.neff_x_max);
  ops >> Option('v', "neff-y-min", opts_.neff_y_min, opts_.neff_y_min);
  ops >> Option('V', "neff-y-max", opts_.neff_y_max, opts_.neff_y_max);
  ops >> Option('y', "neff-y-target", opts_.neff_y_target, opts_.neff_y_target);
  ops >> Option('j', "neff-d-min", opts_.neff_d_min, opts_.neff_d_min);
  ops >> Option('J', "neff-d-max", opts_.neff_d_max, opts_.neff_d_max);
  ops >> OptionPresent('C', "central-mod", opts_.central_mod);
  ops >> Option('R', "round", opts_.round, opts_.round);
  ops >> OptionPresent('M', "use-min-round", opts_.use_min_round);
  ops >> Option('r', "seed", opts_.seed, opts_.seed);

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
  fprintf(out_, "  %-30s %s (def=%s)\n", "-d, --dir <dir>",
          "Directory with count profiles (and A3M alignments)", opts_.dir.c_str());
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file with sampled training set");
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
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-g, --singletons [0,1]",
          "Fraction of singletons in training profiles", opts_.singletons);

  fprintf(out_, "  %-30s %s\n", "-D, --context-data <file>",
          "Context data for adding cs-pseudocounts instead of subst. matrix PCs");
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "-x, --pc-admix [0,1]",
          "Pseudocount admix to be added to count profiles", opts_.pc_admix);
  fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
          "Constant in pseudocount calculation for alignments", opts_.pc_ali);


  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-u, --neff-x-min [1,20]",
          "Minimum number of effective sequences in profiles", opts_.neff_x_min);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-U, --neff-x-max [1,20]",
          "Maximum number of effective sequences in profiles", opts_.neff_x_max);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-v, --neff-y-min [1,20]",
          "Minimum number of effective sequences in pseudocounts column", opts_.neff_y_min);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-V, --neff-y-max [1,20]",
          "Maximum number of effective sequences in pseudocounts column", opts_.neff_y_max);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-y, --neff-y-target [1,inf[",
          "Target number of effective sequences in training pseudocounts column", opts_.neff_y_target);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-j, --neff-d-min [1,20]",
          "Minimum difference of Neff between pseudocounts column and profile", opts_.neff_d_min);
  fprintf(out_, "  %-30s %s (def=%.1f)\n", "-J, --neff-d-max [1,20]",
          "Maximum difference of Neff between pseudocounts column and profile", opts_.neff_d_max);
  fprintf(out_, "  %-30s %s (def=false)\n", "-C, --central-mod", "Modify central profile columns");
  fprintf(out_, "  %-30s %s (def=off)\n", "-R, --round [1,inf[",
          "Use profile ID_R.prf of PSI-BLAST round R");
  fprintf(out_, "  %-30s %s (def=off)\n", "-M, --use-min-round",
          "If several profiles meet the filtering condition, use the one of the minimum round");

  fprintf(out_, "  %-30s %s (def=%u)\n", "-r, --seed [0,inf[",
          "Seed for random number generator", opts_.seed);
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
  ProgressBar progress(out_, 70);

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

  // Read profiles which meet the filtering conditions
  ReadProfiles(progress);

  // Add pseudocounts to count profiles in parallel
  if (opts_.pc_admix > 0.0) {
    fputs("Adding pseudocounts to profiles ...\n", out_);
    CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
    const int nprof = profiles_x_.size();
    progress.Init(nprof);
#pragma omp parallel for schedule(static)
    for (int n = 0; n < nprof; ++n) {
      profiles_x_[n].counts = pc_->AddTo(profiles_x_[n], admix);
      Normalize(profiles_x_[n].counts, profiles_x_[n].neff);
      if (profiles_xy_) {
        CountProfile<Abc>& prof = profiles_y_->at(n);
        prof.counts = pc_->AddTo(prof, admix);
        Normalize(prof.counts, prof.neff);
      }
#pragma omp critical (pseudocounts_advance_progress)
      progress.Advance();
    }
    fputs("\n", out_);
  }

  // Initialize data for computing some statistics
  stats_.reserve(profiles_x_.size());
  for (size_t i = 0; i < profiles_x_.size(); ++i)
    stats_.push_back(vector<size_t>(profiles_x_[i].length(), 0));

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

  // Finally we calculate some sampling statistics
  ComputeStatistics(progress);

  if (profiles_xy_) delete profiles_y_;
  fputs("\nDone!\n", out_);

  return 0;
}

template<class Abc>
void CSTrainSetApp<Abc>::ReadProfiles(ProgressBar& progress) {
  Files files;
  fputs("Globbing input directory for profiles ...", out_);
  fflush(out_);
  GetAllFiles(opts_.dir, files, opts_.profile_ext);
  
  // Group count-profile files of different PSI-BLAST rounds
  typedef string GroupKey;
  typedef std::pair<size_t, size_t> GroupValue;
  typedef std::map<GroupKey, GroupValue> Groups;
  Groups groups;
  // Groups hash:
  //  - count-profile refers to a specific round:
  //    basename = ID_X.prf where X = [1-9]
  //    # key = ID
  //    round = X
  //  - count-profile does not refer to a specifc round:
  //    basename = ID.prf
  //    # key = ID.prf
  //    round = 0
  //  # value = (rounds encoded as bit string, highest round)
  for (Files::iterator it = files.begin(); it != files.end(); ++it) {
    string basename = GetBasename(*it);
    string key = basename;
    size_t round = 0;
    string name = GetBasename(*it, false);
    size_t i = name.rfind('_');
    if (i != string::npos) {
      round = atoi(name.substr(i + 1).c_str());
      if (round > 0) key = name.substr(0, i); // count-profiles refers to a specific round
    }
    Groups::iterator it_groups = groups.find(key);
    if (it_groups == groups.end()) 
      groups[key] = GroupValue(1 << round, round);
    else {
      it_groups->second.first |= 1 << round;
      it_groups->second.second = MAX(it_groups->second.second, round);
    }
  }
  fprintf(out_, "\n%zu files globbed\n", groups.size());
  if (groups.size() == 0) exit(0);
  fprintf(out_, "Reading %zu profiles into memory ...\n", groups.size());
  progress.Init(groups.size());

  // Shuffle group keys
  vector<GroupKey> group_keys;
  group_keys.reserve(groups.size());
  for (Groups::iterator it = groups.begin(); it != groups.end(); ++it)
    group_keys.push_back(it->first);
  Ran ran(opts_.seed);
  random_shuffle(group_keys.begin(), group_keys.end(), ran);

  // profiles_xy_ == true <=> profiles_x_ != profiles_y_ 
  profiles_xy_ = opts_.neff_d_max > 0 && opts_.round == 0;
  if (profiles_xy_)
    profiles_y_ = new ProfileSet();
  else
    profiles_y_ = &profiles_x_;

  for (vector<GroupKey>::iterator it = group_keys.begin(); it != group_keys.end(); ++it) {
    GroupKey key = *it;
    GroupValue value = groups[*it];
    ProfileSet profiles_x;
    ProfileSet profiles_y;
    Files files_x;
    vector<double> profiles_x_neff;
    vector<double> profiles_y_neff;
    double neff_x_min = DBL_MAX;
    double neff_x_max = -DBL_MAX;
    // Read profile of each round
    for (size_t r = 0; r <= value.second; ++r) {
      if (((value.first & (1 << r)) == 0) || 
          (opts_.round > 0 && r != opts_.round)) continue;
      char s[1000];
      if (value.second == 0)
        sprintf(s, "%s/%s", opts_.dir.c_str(), key.c_str());
      else
        sprintf(s, "%s/%s_%zu.%s", opts_.dir.c_str(), key.c_str(), r, opts_.profile_ext.c_str());
      string filename = s;
      FILE* fin = fopen(filename.c_str(), "r");
      if (!fin) throw Exception("Unable to open file '%s'!", filename.c_str());
      CountProfile<Abc> cp(fin);
      fclose(fin);

      // Append profile to profiles_x and profiles_y depending on the filtering contitions
      if (cp.length() < opts_.wlen) continue;        
      double neff = Neff(cp);
      if (!profiles_xy_) {
        // No specific profiles for pseudocounts column
        if (neff >= opts_.neff_x_min && neff <= opts_.neff_x_max &&
            neff >= opts_.neff_y_min && neff <= opts_.neff_y_max) {
          profiles_x.push_back(cp);
          files_x.push_back(filename);
          if (opts_.use_min_round) break;
        }
      } else {
        // Specific profiles for pseudocounts column
        if (neff >= opts_.neff_x_min && neff <= opts_.neff_x_max && 
            (!opts_.use_min_round || profiles_x.size() == 0)) {
          profiles_x.push_back(cp);
          profiles_x_neff.push_back(neff);
          files_x.push_back(filename);
          neff_x_min = MIN(neff_x_min, neff);
          neff_x_max = MAX(neff_x_max, neff);
        } 
        if (neff >= opts_.neff_y_min && neff <= opts_.neff_y_max &&
            (!opts_.use_min_round || profiles_y.size() == 0)) {
          profiles_y.push_back(cp);
          profiles_y_neff.push_back(neff);
        }
        if (opts_.use_min_round && profiles_x.size() > 0 && profiles_y.size() > 0) break;
      }
    }
    progress.Advance();
    if (profiles_x.size() == 0) continue;

    if (!profiles_xy_) {
      // Sample profile from profiles_x
      size_t ix = ran(profiles_x.size());
      profiles_x_.push_back(profiles_x[ix]);
      files_x_.push_back(files_x[ix]);
    } else {
      // Filter profiles_y to meet the distance condition neff_d_(min|max)
      if (profiles_y.size() == 0) continue;
      double neff_y_min = neff_x_min + opts_.neff_d_min;
      double neff_y_max = neff_x_max + opts_.neff_d_max;
      ProfileSet profiles_tmp;
      vector<double> profiles_tmp_neff;
      for (size_t i = 0; i < profiles_y.size(); ++i) {
        if (profiles_y_neff[i] >= neff_y_min && profiles_y_neff[i] <= neff_y_max) {
          profiles_tmp.push_back(profiles_y[i]);
          profiles_tmp_neff.push_back(profiles_y_neff[i]);
        }
      }
      if (profiles_tmp.size() == 0) continue;
      // Sample profile y the for pseudocounts column
      profiles_y = profiles_tmp;
      profiles_y_neff = profiles_tmp_neff;
      size_t iy = ran(profiles_y.size());
      CountProfile<Abc>& profile_y = profiles_y[iy];

      // Filter profiles_x to meet the distance condition neff_d_(min|max)
      double neff_x_min = profiles_y_neff[iy] - opts_.neff_d_max;
      double neff_x_max = profiles_y_neff[iy] - opts_.neff_d_min;
      Files files_tmp;
      profiles_tmp.clear();
      profiles_tmp_neff.clear();
      for (size_t i = 0; i < profiles_x.size(); ++i) {
        if (profiles_x_neff[i] >= neff_x_min && profiles_x_neff[i] <= neff_x_max) {
          profiles_tmp.push_back(profiles_x[i]);
          profiles_tmp_neff.push_back(profiles_x_neff[i]);
          files_tmp.push_back(files_x[i]);
        }
      }
      if (profiles_tmp.size() == 0) continue;
      profiles_x = profiles_tmp;
      profiles_x_neff = profiles_tmp_neff;
      files_x = files_tmp;
      // Sample profile x
      size_t ix = ran(profiles_x.size());
      CountProfile<Abc>& profile_x = profiles_x[ix];

      if (profile_x.length() != profile_y.length())
        throw Exception("Profile sizes differ: '%s'!", files_x[ix].c_str());
      assert(Neff(profile_x) >= opts_.neff_x_min && Neff(profile_x) <= opts_.neff_x_max);
      assert(Neff(profile_y) >= opts_.neff_y_min && Neff(profile_y) <= opts_.neff_y_max);
      assert(Neff(profile_y) - Neff(profile_x) >= opts_.neff_d_min && 
           Neff(profile_y) - Neff(profile_x) <= opts_.neff_d_max);

      profiles_x_.push_back(profile_x);
      profiles_y_->push_back(profile_y);
      files_x_.push_back(files_x[ix]);
    }
  }
  assert(profiles_x_.size() == profiles_y_->size());
  assert(profiles_x_.size() == files_x_.size());
  fprintf(out_, "\n%zu profiles passed filter\n", profiles_x_.size());
  if (profiles_x_.size() == 0) exit(0);
}

template<class Abc>
void CSTrainSetApp<Abc>::ComputeStatistics(ProgressBar& progress) {
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
}

template<class Abc>
void CSTrainSetApp<Abc>::SampleFullLengthProfiles(ProfileSet& samples) {
  const size_t nsamples = MIN(profiles_x_.size(), opts_.nsamples);
  fprintf(out_, "Sampling %zu profiles ...\n", nsamples);
  ProgressBar progress(out_, 70, nsamples);

  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    if (samples.size() < nsamples) {
      samples.push_back(profiles_x_[n]);
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
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    const CountProfile<Abc>& cp = profiles_x_[n];
    const size_t ncol =  cp.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (cp.neff[i + center] >= opts_.neff_x_min &&
          cp.neff[i + center] <= opts_.neff_x_max) unmasked++;
    }
    nall += MIN(opts_.max_win, unmasked);
  }
  const size_t nsamples = MIN(opts_.nsamples, nall);

  fprintf(out_, "Sampling %zu profiles with W=%zu out of %zu windows ...\n",
          nsamples, opts_.wlen, nall);

  // Iterate over input data and build counts profiles
  ProgressBar progress(out_, 70, nsamples);
#pragma omp parallel for schedule(static)
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    if (samples.size() < nsamples) {
      Ran ran(opts_.seed + n);
      const CountProfile<Abc>& cp = profiles_x_[n];
      LOG(ERROR) << "sampling " << files_x_[n];

      vector<int> shuffle;
      for (size_t i = 0; i <= cp.counts.length() - opts_.wlen; ++i)
        if (cp.neff[i + center] >= opts_.neff_x_min &&
            cp.neff[i + center] <= opts_.neff_x_max) shuffle.push_back(i);
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
inline bool CSTrainSetApp<Abc>::IsValidPos(double neff_x, double neff_y) const {
  return neff_x >= opts_.neff_x_min && neff_x <= opts_.neff_x_max &&
         neff_y >= opts_.neff_y_min && neff_y <= opts_.neff_y_max &&
         neff_y - neff_x >= opts_.neff_d_min &&
         neff_y - neff_x <= opts_.neff_d_max;
}

template<class Abc>
void CSTrainSetApp<Abc>::SampleTrainingSeqs(TrainSeqs& samples, size_t nsamples_) {
  const size_t center = (opts_.wlen - 1) / 2;
  Vector<double> neff(profiles_x_.size(), 0.0);
  double nall = 0.0;  // total pool size

  // Precompute Neff and total pool size taking masking into account
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    const CountProfile<Abc>& cp_x = profiles_x_[n];
    const size_t ncol = cp_x.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (IsValidPos(cp_x.neff[i + center], profiles_y_->at(n).neff[i + center])) {
        unmasked++;
        neff[n] += opts_.sampling_mode == 2 ? cp_x.neff[i + center] : 1.0;
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
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    if (samples.size() < nsamples) {
      Ran ran(opts_.seed + n);
      const CountProfile<Abc>& cp_x = profiles_x_[n];
      size_t nmatch = cp_x.counts.length();

      // Filter out columns that don't suffice Neff-filter, same as above
      size_t ncol = nmatch - opts_.wlen + 1;
      Vector<bool> masked(ncol, false);
      int unmasked = 0;
      for (size_t i = 0; i < ncol; ++i) {
        if (IsValidPos(cp_x.neff[i + center], profiles_y_->at(n).neff[i + center]))
          unmasked++;
        else
          masked[i] = true;
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
      string file = opts_.dir + kDirSep + GetBasename(files_x_[n], false) + "." + opts_.ali_ext;
      FILE* fp = fopen(file.c_str(), "r");
      if (!fp) throw Exception("Can't open alignment file '%s'!", file.c_str());
      Alignment<Abc> ali(fp, A3M_ALIGNMENT);
      fclose(fp);

      // Check if number of match columns in profile and alignment are correct
      if (nmatch != ali.nmatch())
        throw Exception("Number of matchcols in ali '%s' should be %zu but is %zu!",
                        GetBasename(file).c_str(), cp_x.counts.length(), nmatch);

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
        double neff_center = opts_.sampling_mode == 2 ? cp_x.neff[m + center] : 1.0;
        double s = 1.1 * todo * (neff_center / left);
        LOG(ERROR) << strprintf("todo=%5.2f  neff[%2zu]=%4.1f  left=%9.1f  s=%4.2f",
                                todo, m, cp_x.neff[m + center], left, s);
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
              // Cut training pseudocounts column from profiles_y_
              CountProfile<Abc> cpw_y = CountProfile<Abc>(profiles_y_->at(n), m, opts_.wlen);
              if (opts_.neff_y_target != 0.0 && cpw_y.neff[center] < opts_.neff_y_target) {               
                Profile<Abc> pw_y = cpw_y.counts;
                Normalize(pw_y, 1.0);
                // Estimated target Neff in count profile = Neff / Neff(cp_i) * Neff(p)
                CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
                pc_->SetTargetNeff(opts_.neff_y_target / cpw_y.neff[center] * Neff(pw_y));
                cpw_y.counts = pc_->AddTo(cpw_y, admix);
                Assign(cpw_y.neff, opts_.neff_y_target);
                Normalize(cpw_y.counts, cpw_y.neff);
              }
              ProfileColumn<Abc> y(cpw_y.counts[center]);
              if (samples.size() < nsamples) {
#pragma omp critical (add_sample)
                {
                  samples.push_back(TrainingSequence<Abc>(seq, y));
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
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    const CountProfile<Abc>& cp_x = profiles_x_[n];
    const size_t ncol =  cp_x.counts.length() - opts_.wlen + 1;
    size_t unmasked = 0;
    for (size_t i = 0; i < ncol; ++i) {
      if (IsValidPos(cp_x.neff[i + center], profiles_y_->at(n).neff[i + center]))
        unmasked++;
    }
    nall += MIN(opts_.max_win, unmasked);
  }
  const size_t nsamples = MIN(opts_.nsamples - nsingletons, nall);

  fprintf(out_, "Sampling %zu training profiles with W=%zu out of %zu windows ...\n",
          nsamples, opts_.wlen, nall);

  // Iterate over input data and build training profiles
  ProgressBar progress(out_, 70, nsamples);
#pragma omp parallel for schedule(static)
  for (size_t n = 0; n < profiles_x_.size(); ++n) {
    if (samples.size() < nsamples + nsingletons) {
      Ran ran(opts_.seed + n + 1);
      const CountProfile<Abc>& cp_x = profiles_x_[n];
      LOG(ERROR) << "sampling " << files_x_[n];

      vector<int> shuffle;
      for (size_t i = 0; i <= cp_x.counts.length() - opts_.wlen; ++i) {
        if (IsValidPos(cp_x.neff[i + center], profiles_y_->at(n).neff[i + center]))
          shuffle.push_back(i);
      }
      random_shuffle(shuffle.begin(), shuffle.end(), ran);

      // Proportional number of training windows we need to sample
      size_t w = MIN(opts_.max_win, shuffle.size());
      size_t s = ceil(nsamples * (static_cast<double>(w) / nall));
      LOG(ERROR) << strprintf("n=%zu nsamples=%zu  nall=%zu  len=%zu  s=%zu",
                              n, nsamples, nall, cp_x.counts.length(), s);
      s = MIN(s, shuffle.size());

      // Copy profile windows into 'samples' vector
      for (size_t i = 0; i < s; ++i) {
        // Cut the profile window
        CountProfile<Abc> cpw_x(cp_x, shuffle[i], opts_.wlen);
        if (opts_.central_mod) {
          double* central = cpw_x.counts[center];
          size_t max = 0;
          for (size_t a = 1; a < Abc::kSize; ++a)
            if (central[a] > central[max]) max = a;
          for (size_t a = 0; a < Abc::kSize; ++a)
            central[a] = a == max ? cpw_x.neff[center] : 0.0;
        }
        // Cut the profile column
        CountProfile<Abc> cpw_y = CountProfile<Abc>(profiles_y_->at(n), shuffle[i], opts_.wlen);
        if (opts_.neff_y_target != 0.0 && cpw_y.neff[center] < opts_.neff_y_target) {               
          Profile<Abc> pw_y = cpw_y.counts;
          Normalize(pw_y, 1.0);
          // Estimated target Neff in count profile = Neff / Neff(M_i) * Neff(CP_i)
          CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
          pc_->SetTargetNeff(opts_.neff_y_target / cpw_y.neff[center] * Neff(pw_y));
          pw_y = pc_->AddTo(cpw_y, admix);
          cpw_y.counts = pw_y;
          Assign(cpw_y.neff, opts_.neff_y_target);
          Normalize(cpw_y.counts, cpw_y.neff);
        }
        ProfileColumn<Abc> y(cpw_y.counts[center]);
        // Add the training profile
#pragma omp critical (add_sample)
        if (samples.size() < nsamples + nsingletons) {
          Normalize(cpw_x.counts, cpw_x.neff);
          samples.push_back(TrainingProfile<Abc>(cpw_x, y));
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

