// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "context_library-inl.h"
#include "count_profile-inl.h"
#include "emission.h"
#include "progress_bar.h"
#include "sequence-inl.h"

using namespace GetOpt;
using std::string;
using std::vector;
using std::map;

namespace cs {

struct CSMakematAppOptions {
  CSMakematAppOptions() { Init(); }
  virtual ~CSMakematAppOptions() {}

  // Set csbuild default parameters
  void Init() {
    conf = 8;
    rmin = 0;
    rmax = 3;
    weight_center = 1.6;
    weight_decay  = 0.85;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (queryfile.empty()) throw Exception("No query file provided!");
    if (dbfile.empty()) throw Exception("No database file provided!");
    if (libfile.empty()) throw Exception("No context library provided!");
  }

  string infile;     // input file with aligned pairs and confidence values
  string outfile;    // output file for subst. matrix target probabilities
  string queryfile;  // input file with concatenated count profiles as queries
  string dbfile;     // input file with FASTA formatted abstract state sequences
  string libfile;    // input file with context profile library
  int conf;          // minimum posterior value for inclusion in counts
  int rmin;          // include query profiles starting at this HHblits round
  int rmax;          // include query profiles up to this HHblits round
  double weight_center;  // weight of central column in multinomial emission
  double weight_decay;   // exponential decay of window weights
};  // CSMakematAppOptions


class CSMakematApp : public Application {
 protected:
  typedef map<string, CountProfile<AA> > ProfilesMap;
  typedef map<string, Sequence<AS62> > SequencesMap;

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

  // Parameter wrapper
  CSMakematAppOptions opts_;
};  // class CSMakematApp


void CSMakematApp::ParseOptions(GetOpt_pp& ops) {
  string rounds = "";
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option('q', "queryfile", opts_.queryfile, opts_.queryfile);
  ops >> Option('d', "dbfile", opts_.dbfile, opts_.dbfile);
  ops >> Option('c', "conf", opts_.conf, opts_.conf);
  ops >> Option('r', "rounds", rounds, rounds);
  ops >> Option('D', "context-data", opts_.libfile, opts_.libfile);

  // Parse rmin,rmax option
  if (!rounds.empty()) {
    if (rounds.find(',') != string::npos) {  // rmin,rmax given
      vector<string> rmin_rmax;
      Split(rounds, ',', rmin_rmax);
      opts_.rmin = atoi(rmin_rmax[0].c_str());
      opts_.rmax = atoi(rmin_rmax[1].c_str());
    } else {  // rmin=rmax given
      opts_.rmin = atoi(rounds.c_str());
      opts_.rmax = opts_.rmin;
    }
  }

  opts_.Validate();

  if (opts_.outfile.empty())
    opts_.outfile = GetBasename(opts_.infile, false) + "mat";
}

void CSMakematApp::PrintBanner() const {
  fputs("Construct an abstract state substit. matrix from aligned pairs.\n", out_);
}

void CSMakematApp::PrintUsage() const {
  fputs("Usage: csmakemat -i <pairs> -q <queries> -d <db> -D <lib> [options]\n",
        out_);
}

void CSMakematApp::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input file with aligned pairs and confidences");
  fprintf(out_, "  %-30s %s\n", "-q, --queryfile <file>",
          "Input file with query profiles");
  fprintf(out_, "  %-30s %s\n", "-d, --dbfile <file>",
          "Input file with abstract state database sequences");
  fprintf(out_, "  %-30s %s\n", "-D, --context-data <file>",
          "Input file with context library");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Output file for substitution matrix");
  fprintf(out_, "  %-30s %s (def=%d)\n", "-c, --conf [0:9]",
          "Minimum confidence value of aligned pairs", opts_.conf);
  fprintf(out_, "  %-30s %s (def=%d,%d)\n", "-r, --rounds [0:9],[0:9]",
          "Range of HHblits rounds to be included", opts_.rmin, opts_.rmax);
}

int CSMakematApp::Run() {
  // First we read the context library
  FILE* fin = fopen(opts_.libfile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.libfile.c_str());
  ContextLibrary<AA> lib(fin);
  TransformToLog(lib);
  fclose(fin);
  Emission<AA> emission(lib.wlen(), opts_.weight_center, opts_.weight_decay);

  // Read all querry profiles and store in map
  ProfilesMap query_profiles;
  fputs("Reading abstract state query profiles ...", out_);
  fflush(out_);
  fin = fopen(opts_.queryfile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.queryfile.c_str());
  size_t num_queries = 0;  // number of unique queries without HHblits rounds
  while (!feof(fin)) {
    // Read next count profile
    CountProfile<AA> cp(fin);
    query_profiles.insert(std::make_pair(cp.name, cp));
    if (cp.name.find('_') == string::npos) num_queries++;
    // Check for EOF
    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
  fclose(fin);
  fprintf(out_, " %zu profiles read \n", query_profiles.size());

  // Read all abstract state sequences and store in map
  SequencesMap database_seqs;
  fputs("Reading abstract state database sequences ...", out_);
  fflush(out_);
  fin = fopen(opts_.dbfile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.dbfile.c_str());
  vector<string> tokens;
  while (!feof(fin)) {
    // Read next sequence
    Sequence<AS62> seq(fin);
    Split(seq.header(), '|', tokens);
    database_seqs.insert(std::make_pair(tokens[1], seq));
    tokens.clear();

    // Check for EOF
    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
  fclose(fin);
  fprintf(out_, " %zu sequences read \n", database_seqs.size());

  // Parse file with aligned pairs from HHsearch results and count substitutions
  vector<const CountProfile<AA>* > profiles;  // pointers to profiles of each round
  const Sequence<AS62>* seq = NULL;           // pointer to abstract state sequence
  vector<Matrix<double> > probs;              // posterior matrices of each round
  Matrix<double> counts(lib.size(), AS62::kSize, 0.0);  // substitution counts
  string qname, tname;
  ProgressBar progress(out_, 70, num_queries);
  const size_t kBuffSize = 16 * KB;
  char buffer[kBuffSize];

  fputs("Counting context to abstract state substitutions ...\n", out_);
  fin = fopen(opts_.infile.c_str(), "r");
  if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());

  while (fgetline(buffer, kBuffSize, fin)) {
    if (!strscn(buffer)) continue;
    Split(buffer, '\t', tokens);

    if (tokens.size() == 2) {  // is this line a new hit block?
      if (tokens[0] != qname) {  // new query or same query with new hit?
        // Fetch sequence based abstract state profile
        profiles.clear();
        if (opts_.rmin == 0) {
          qname = tokens[0];
          ProfilesMap::iterator pi = query_profiles.find(qname);
          if(pi != query_profiles.end()) profiles.push_back(&(pi->second));
          else throw Exception("Cannot find query profile '%s'!", qname.c_str());
        }

        // Fetch abstract state profiles of HHblits rounds _1, _2, etc.
        for (int r = MAX(1, opts_.rmin); r <= opts_.rmax; ++r) {
          qname = tokens[0] + strprintf("_%d", r);
          ProfilesMap::iterator pi = query_profiles.find(qname);
          if(pi != query_profiles.end()) profiles.push_back(&(pi->second));
          else break;
        }
        qname = tokens[0];

        // For each round calculate posterior probabilities given its query profile
        probs.clear();
        for (size_t r = 0; r < profiles.size(); ++r) {
          Matrix<double> pp(profiles[r]->length(), lib.size(), 0.0);
          for (size_t i = 0; i < profiles[r]->length(); ++i)
            CalculatePosteriorProbs(lib, emission, *profiles[r], i, &pp[i][0]);
          probs.push_back(pp);
        }

        progress.Advance(); // advance progress
      }

      // Fetch database side abstract state sequence of hit
      tname = tokens[1];
      SequencesMap::iterator si = database_seqs.find(tname);
      if(si != database_seqs.end()) seq = &(si->second);
      else throw Exception("Cannot find database sequence '%s'!", tname.c_str());
      LOG(INFO) << seq->header();

    } else if (qname != tname && atoi(tokens[2].c_str()) >= opts_.conf) {
      // Line contains an aligned residue pair with sufficient confidence value
      size_t i = atoi(tokens[0].c_str()) - 1;
      size_t j = atoi(tokens[1].c_str()) - 1;
      assert(j < seq->length());
      assert(seq != NULL);

      for (size_t r = 0; r < profiles.size(); ++r) {
        assert(i < profiles[r]->length());
        assert(profiles[r] != NULL);
        size_t l = (*seq)[j];
        for (size_t k = 0; k < lib.size(); ++k)
          counts[k][l] += probs[r][i][k];
      }
    }

    tokens.clear();
  }
  progress.Complete();
  fclose(fin);
  fputs("\n", out_);

  // Calculate sum, min, and max over all cells
  double sum = 0.0, min = FLT_MAX, max = 0.0;
  for (size_t k = 0; k < lib.size(); ++k) {
    for (size_t l = 0; l < AS62::kSize; ++l) {
      sum += counts[k][l];
      min = MIN(min, counts[k][l]);
      max = MAX(max, counts[k][l]);
    }
  }
  fputs("Count statistics:\n", out_);
  fprintf(out_, "Sum  = %.2f\n", sum);
  fprintf(out_, "Min  = %.2f\n", min);
  fprintf(out_, "Max  = %.2f\n", max);
  fprintf(out_, "Mean = %.2f\n", sum / (lib.size() * AS62::kSize));

  // Calculate target frequencies for abstract state substitution matrix
  Matrix<double> q(lib.size(), AS62::kSize, 0.0);
  for (size_t k = 0; k < lib.size(); ++k) {
    for (size_t l = 0; l < AS62::kSize; ++l) {
      q[k][l] = counts[k][l] / sum;
    }
  }

  // Write matrix to outfile
  FILE* fout = fopen(opts_.outfile.c_str(), "w");
  if (!fout) throw Exception("Cannot write to file '%s'!", opts_.outfile.c_str());
  for (size_t k = 0; k < AS62::kSize; ++k)
    fprintf(fout, "%c%7c ", k== 0 ? '#' : ' ', AS62::kIntToChar[k]);
  fputs("\n", fout);
  for (size_t k = 0; k < lib.size(); ++k) {
    for (size_t l = 0; l < AS62::kSize; ++l) {
      fprintf(fout, "%8.2g ", q[k][l]);
    }
    fputs("\n", fout);
  }
  fclose(fout);

  fprintf(out_, "Wrote target freqs of abstract state substitution matrix to %s\n",
          opts_.outfile.c_str());
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  return cs::CSMakematApp().main(argc, argv, stdout, "csmakemat");
}
