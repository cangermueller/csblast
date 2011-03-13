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

using namespace GetOpt;
using std::string;
using std::vector;
using std::count;

namespace cs {

struct CSTrainSetAppOptions {
    CSTrainSetAppOptions() { Init(); }

    void Init() {
        dir           = "./";
        nsamples      = 1000000;
        wlen          = 0;
        seed          = 0;
        pc_admix      = 0.01;
        pc_ali        = 12.0;
        min_neff      = 0;
        max_win       = 1000;
        sampling_mode = 0;
        weight_center = 1.6;
        weight_decay  = 0.85;
        pc_engine     = "auto";
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (sampling_mode > 0 && wlen == 0) throw Exception("No window length given!");
        if (outfile.empty()) throw Exception("No output file provided!");
        if (wlen > 0 && !(wlen & 1)) throw Exception("Window length must be odd!");
    }

    string dir;            // directory with input data to build trainset from
    string ext;            // file extension of input files with counts
    string outfile;        // created training set
    string modelfile;      // input file with context profile library or HMM
    string pc_engine;      // pseudocount engine
    size_t nsamples;       // size of training set to create
    size_t wlen;           // length of context window (zero means full-length)
    double pc_admix;       // minimal pseudocount admix for count profiles
    double pc_ali;         // constant in pseudocount calculation for alignments
    double min_neff;       // minimum Neff of count profiles
    size_t max_win;        // maximal number of windows per protein
    int sampling_mode;     // sample sequence-column pairs from alignments
    double weight_center;  // weight of central column in multinomial emission
    double weight_decay;   // exponential decay of window weights
    unsigned int seed;     // seed
};  // CSTrainSetAppOptions


template<class Abc>
class CSTrainSetApp : public Application {
  protected:
    typedef vector<CountProfile<Abc> > ProfileSet;
    typedef vector<TrainingSequence<Abc> > TrainSeqs;
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
    void SampleTrainingSeqs(TrainSeqs& samples);
    // Initializes substitution matrix (specialized by alphabet type).
    void InitSubstitutionMatrix();

    CSTrainSetAppOptions opts_;
    ProfileSet profiles_;
    Files files_;
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;
    scoped_ptr<ContextLibrary<Abc> > lib_;
    scoped_ptr<Crf<Abc> > crf_;
    scoped_ptr<Pseudocounts<Abc> > pc_;
};  // CSTrainSetApp


template<class Abc>
void CSTrainSetApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
    ops >> Option('d', "dir", opts_.dir, opts_.dir);
    ops >> Option('N', "size", opts_.nsamples, opts_.nsamples);
    ops >> Option('W', "wlen", opts_.wlen, opts_.wlen);
    ops >> Option('r', "seed", opts_.seed, opts_.seed);
    ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
    ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
    ops >> Option('n', "neff", opts_.min_neff, opts_.min_neff);
    ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
    ops >> Option('s', "sampling-mode", opts_.sampling_mode, opts_.sampling_mode);

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
    fprintf(out_, "  %-30s %s (def=%zu)\n", "-N, --size [0,inf[",
            "Size of training set to sample", opts_.nsamples);
    fprintf(out_, "  %-30s %s\n", "-W, --wlen [0,inf[",
            "Sample profile windows of length W instead of full-length");
    fprintf(out_, "  %-30s %s (def=%d)\n", "-s, --sampling-mode 0|1|2",
            "Sampling mode", opts_.sampling_mode);
    fprintf(out_, "  %-30s %s\n", "", "0: sample count profiles");
    fprintf(out_, "  %-30s %s\n", "", "1: sample training seqs from alignment query");
    fprintf(out_, "  %-30s %s\n", "", "2: sample training seqs from full alignment");
    fprintf(out_, "  %-30s %s\n", "-D, --context-data <file>",
            "Context data for adding cs-pseudocounts instead of subst. matrix PCs");
    fprintf(out_, "  %-30s %s (def=%.2f)\n", "-x, --pc-admix [0,1]",
            "Pseudocount admix to be added to count profiles", opts_.pc_admix);
    fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
            "Constant in pseudocount calculation for alignments", opts_.pc_ali);
    fprintf(out_, "  %-30s %s (def=%.1f)\n", "-n, --neff [1,20]",
            "Minimum number of effective sequences in profiles", opts_.min_neff);
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

    // Glob directory for count-profile files with .prf extension
    fputs("Globbing input directory for profiles ...", out_);
    fflush(out_);
    GetAllFiles(opts_.dir, files_, "prf");
    random_shuffle(files_.begin(), files_.end(), ran);
    fprintf(out_, " %zu files globbed\n", files_.size());

    // Read all profiles and calcualte total number of residues
    fprintf(out_, "Reading %zu profiles into memory ...\n", files_.size());
    ProgressBar progress(out_, 70, files_.size());
    Files::iterator it = files_.begin();
    while (it != files_.end()) {
        // Read next count profile
        string filename = opts_.dir + kDirSep + *it;
        LOG(ERROR) << *it;
        FILE* fp = fopen(filename.c_str(), "r");
        if (!fp) throw Exception("Unable to open file '%s'!", filename.c_str());
        CountProfile<Abc> prof(fp);
        fclose(fp);
        progress.Advance();

        // Filter out count profiles with insufficient Neff
        if (Neff(prof) < opts_.min_neff || prof.counts.length() < opts_.wlen) {
            it = files_.erase(it);
            continue;
        }

        // Add profile to pool
        profiles_.push_back(prof);
        it++;
    }
    assert_eq(profiles_.size(), files_.size());
    fprintf(out_, "\n%zu profiles passed filter for minimum Neff=%.2f\n",
            profiles_.size(), opts_.min_neff);

    // Add pseudocounts to count profiles in parallel
    fputs("Adding pseudocounts to profiles ...\n", out_);
    CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
    const int nprof = profiles_.size();
    progress.Init(nprof);
#pragma omp parallel for schedule(static)
    for (int n = 0; n < nprof; ++n) {
        profiles_[n].counts = pc_->AddTo(profiles_[n], admix);
        Normalize(profiles_[n].counts, profiles_[n].neff);
#pragma omp critical (pseudocounts_advance_progress)
        progress.Advance();
    }
    fputs("\n", out_);

    // Now we can call the actual sampling method
    if (opts_.sampling_mode > 0) {
        TrainSeqs samples;
        samples.reserve(opts_.nsamples);
        SampleTrainingSeqs(samples);

        FILE* fout = fopen(opts_.outfile.c_str(), "w");
        if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
        fprintf(out_, "Writing %zu training sequences to %s ...", samples.size(),
                opts_.outfile.c_str());
        fflush(out_);
        WriteAll(samples, fout);
        fclose(fout);

    } else if (opts_.wlen > 0) {
        ProfileSet samples;
        samples.reserve(opts_.nsamples);
        SampleProfileWindows(samples);

        FILE* fout = fopen(opts_.outfile.c_str(), "w");
        if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
        fprintf(out_, "Writing %zu profile windows to %s ...", samples.size(),
                opts_.outfile.c_str());
        fflush(out_);
        WriteAll(samples, fout);
        fclose(fout);

    } else {
        ProfileSet samples;
        samples.reserve(opts_.nsamples);
        SampleFullLengthProfiles(samples);

        FILE* fout = fopen(opts_.outfile.c_str(), "w");
        if (!fout) throw Exception("Can't open outfile '%s'!", opts_.outfile.c_str());
        fprintf(out_, "Writing %zu full-length profiles to %s ...", samples.size(),
                opts_.outfile.c_str());
        fflush(out_);
        WriteAll(samples, fout);
        fclose(fout);
    }
    fputs("\n", out_);

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

    // Iterate over input data and build counts profiles either by full-length
    // conversion or by sampling of context windows
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
            if (s < shuffle.size()) shuffle.erase(shuffle.begin() + s, shuffle.end());
            // Copy profile windows into 'samples' vector
            for (size_t i = 0; i < shuffle.size(); ++i) {
                CountProfile<Abc> p(cp, shuffle[i], opts_.wlen);
#pragma omp critical (add_sample)
                if (samples.size() < nsamples) {
                    samples.push_back(p);
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
void CSTrainSetApp<Abc>::SampleTrainingSeqs(TrainSeqs& samples) {
    const size_t center = (opts_.wlen - 1) / 2;
    Vector<double> neff(profiles_.size(), 0.0);
    double nall = 0.0;  // total pool size

    // Precompute Neff and total pool size taking masking into account
    for (size_t n = 0; n < profiles_.size(); ++n) {
        const CountProfile<Abc>& cp = profiles_[n];
        const size_t ncol = cp.counts.length() - opts_.wlen + 1;
        size_t unmasked = 0;
        for (size_t i = 0; i < ncol; ++i) {
            if (cp.neff[i + center] >= opts_.min_neff) {
                unmasked++;
                neff[n] += opts_.sampling_mode == 1 ? 1.0 : cp.neff[i + center];
            }
        }
        neff[n] /= unmasked;  // average Neff in columns that passed Neff-filter
        nall += MIN(opts_.max_win, unmasked) * neff[n];
    }
    const size_t nsamples = MIN(opts_.nsamples, floor(nall));

    fprintf(out_, "Sampling %zu seq-col pairs out of %.0f virtual windows ...\n",
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

            // Read-in the underlying A3M alignment of profile 'n'
            string file = opts_.dir + kDirSep + GetBasename(files_[n], false) + "a3m";
            FILE* fp = fopen(file.c_str(), "r");
            if (!fp) throw Exception("Can't open alignment file '%s'!", file.c_str());
            Alignment<Abc> ali(fp, A3M_ALIGNMENT);
            fclose(fp);
            size_t nmatch = ali.nmatch();

            // Check if number of match columns in profile and alignment are correct
            if (nmatch != cp.counts.length())
                throw Exception("Number of matchcols in ali '%s' should be %zu but is %zu!",
                                GetBasename(file).c_str(), cp.counts.length(), nmatch);

            // Filter out columns that don't suffice Neff-filter, same as above
            size_t ncol = nmatch - opts_.wlen + 1;
            Vector<bool> masked(ncol, false);
            int unmasked = 0;
            for (size_t i = 0; i < ncol; ++i) {
                if (cp.neff[i + center] < opts_.min_neff) masked[i] = true;
                else unmasked++;
            }

            // Determine workload for this profile
            double nprf = MIN(opts_.max_win, static_cast<size_t>(unmasked)) * neff[n];
            double todo = MIN(1.0, nprf / MAX(0.0, nall)) * (nsamples - samples.size());
            double left = unmasked * neff[n];

            LOG(ERROR) <<
                strprintf("n=%4zu neff=%4.1f L=%4zu nall=%9.1f nprf=%9.1f todo=%5.1f",
                          n, neff[n], nmatch, nall, nprf, todo);

#pragma omp atomic
            nall -= nprf;

            // Try out sequences within each column in random order
            vector<int> shuffle;
            if (opts_.sampling_mode == 1) {
                shuffle.push_back(0);
            } else {
                for (size_t k = 0; k < ali.nseqs(); ++k) shuffle.push_back(k);
            }

            while (todo > 0.0 && samples.size() < nsamples && unmasked > 0) {
                // Pick a random column as window start position
                size_t m = ran(ncol);
                if (masked[m]) continue;
                // Determine sampling work load for column 'm' (10% more as safe-guard)
                double neff_center = opts_.sampling_mode == 1 ? 1.0 : cp.neff[m + center];
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
                            ProfileColumn<Abc> col(&cp.counts[m + center][0]);
#pragma omp critical (add_sample)
                            {
                                samples.push_back(TrainingSequence<Abc>(seq, col));
                                progress.Advance();
                            }
                            s -= 1;     // reduce work load for this column by onea
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

    fputs("\nShuffling sampled sequence-column pairs ...", out_);
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
