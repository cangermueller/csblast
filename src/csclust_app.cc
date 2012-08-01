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
// #include <clocale>
#include "application.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "em_clustering.h"
#include "getopt_pp.h"
#include "matrix_pseudocounts-inl.h"
#include "pdf_writer-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"
#include "tamura_nei_matrix.h"

using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct CSClustAppOptions {
    CSClustAppOptions() { Init(); }

    void Init() {
        toll        = 1e-4;
        nprofs      = 0;
        blosum_type = "BLOSUM62";
        gauss_init  = 0.0;
        seed        = 0;
        pc_init     = 0.5;
        max_iters   = 100;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (nprofs == 0 && libfile.empty())
            throw Exception("No value for number of profiles provided!");
        if (infile.empty())
            throw Exception("No input file with training data provided!");
    }

    string infile;          // infile with training set of count profiles
    string outfile;         // output file with clustered context library
    string pdffile;         // output file for visualized context library
    string libfile;         // library input file for restarting
    size_t nprofs;          // number of context profiles
    double toll;            // LL change for convergence
    int max_iters;          // maximal number of EM iterations to perform
    double gauss_init;      // sigma for gauus initialization
    double pc_init;         // pseudocount admix for sampling initialization
    string blosum_type;     // BLOSUM matrix for pseudocount generation.
    unsigned int seed;      // seed
    EMClusteringParams em;  // wrapper for EM clustering parameters
};  // CSClustAppOptions


template<class Abc>
class CSClustApp : public Application {
  private:
    typedef vector<CountProfile<Abc> > TrainingSet;

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
    // Reads training data from infile.
    void ReadTrainingData();
    // Initializes HMM from jumpstart file or by seeding from training data.
    void InitLibrary();
    // Initializes substitution matrix (specialized by alphabet type).
    void InitSubstitutionMatrix();

    CSClustAppOptions opts_;
    TrainingSet trainset_;
    scoped_ptr<ContextLibrary<Abc> > lib_;
    scoped_ptr<SubstitutionMatrix<Abc> > sm_;
};  // CSClustApp


template<class Abc>
void CSClustApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "infile", opts_.infile, opts_.infile);
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
    ops >> Option('p', "pdf", opts_.pdffile, opts_.pdffile);
    ops >> Option('K', "profiles", opts_.nprofs, opts_.nprofs);
    ops >> Option('t', "toll", opts_.toll, opts_.toll);
    ops >> Option('j', "jumpstart", opts_.libfile, opts_.libfile);
    ops >> Option('N', "iters", opts_.max_iters, opts_.max_iters);
    ops >> Option('x', "pc-admix", opts_.em.pca, opts_.em.pca);
    ops >> Option('g', "gauss-init", opts_.gauss_init, opts_.gauss_init);
    ops >> Option('r', "seed", opts_.seed, opts_.seed);
    ops >> Option('c', "weight-center", opts_.em.weight_center,
                  opts_.em.weight_center);
    ops >> Option('d', "weight-decay", opts_.em.weight_decay,
                  opts_.em.weight_decay);
    ops >> Option(' ', "pc-init", opts_.pc_init, opts_.pc_init);
    opts_.Validate();

    if (opts_.outfile.empty())
        opts_.outfile = GetDirname(opts_.infile) + GetBasename(opts_.infile, false) + "lib";
}

template<class Abc>
void CSClustApp<Abc>::PrintBanner() const {
    fputs("Cluster a training set of profile windows into a context library.\n",
          out_);
}

template<class Abc>
void CSClustApp<Abc>::PrintUsage() const {
    fputs("Usage: csclust -i <trainsetfile> -o <outfile> -K <profiles> [options]\n", out_);
    fputs("       csclust -i <trainsetfile> -o <outfile> -j <libfile> [options]\n", out_);
}

template<class Abc>
void CSClustApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n",             "-i, --infile <file>", "Infile with training set");
    fprintf(out_, "  %-30s %s\n",             "-o, --outfile <file>", "Output file for library");
    fprintf(out_, "  %-30s %s\n",             "-K, --profiles [0,inf[", "Number of profiles");
    fprintf(out_, "  %-30s %s (def=%3.1g)\n", "-t, --toll [0,inf[", "Log-likelihood change per column for convergence", opts_.toll);
    fprintf(out_, "  %-30s %s (def=%d)\n",    "-N, --iters [0,inf[", "Maximal number of EM iterations", opts_.max_iters);
    fprintf(out_, "  %-30s %s\n",             "-j, --jumpstart <file>", "Jumpstart the clustering with a context library");
    fprintf(out_, "  %-30s %s (def=%.2f)\n",  "-x, --pc-admix [0,inf[", "Pseudocount admix for profile probs in EM clustering", opts_.em.pca);
    fprintf(out_, "  %-30s %s (def=off)\n",   "-g, --gauss-init [0,inf[", "Turn on gaussian profile initialization using given sigma");
    fprintf(out_, "  %-30s %s (def=%u)\n",    "-r, --seed [0,inf[", "Seed for random number generator", opts_.seed);
    fprintf(out_, "  %-30s %s (def=%.2f)\n",  "-p, --pc-init [0,1]", "Pseudocount admix in profile initialization by sampling", opts_.pc_init);
    fprintf(out_, "  %-30s %s (def=%4.2f)\n", "-c, --weight-center [0,inf[", "Weight of central profile column in context window", opts_.em.weight_center);
    fprintf(out_, "  %-30s %s (def=%4.2f)\n", "-d, --weight-decay [0,1]", "Exponential decay of positional window weights", opts_.em.weight_decay);
    fprintf(out_, "  %-30s %s\n",             "    --pdf <file>", "Generate PDF with profile logos of learned context profies.");
}

template<class Abc>
void CSClustApp<Abc>::InitSubstitutionMatrix() {
    BlosumType type = BlosumTypeFromString(opts_.blosum_type);
    sm_.reset(new BlosumMatrix(type));
}

template<>
void CSClustApp<Dna>::InitSubstitutionMatrix() {
    sm_.reset(new TamuraNeiMatrix());
}

template<class Abc>
void CSClustApp<Abc>::ReadTrainingData() {
    FILE* fin;
    fprintf(out_, "Reading training set from %s ...\n",
            GetBasename(opts_.infile).c_str());
    fin = fopen(opts_.infile.c_str(), "r");
    if (!fin)
        throw Exception("Can't read training set from '%s'!", opts_.infile.c_str());
    ReadAll(fin, trainset_);
    fclose(fin);
    fprintf(out_, "%zu profiles read\n", trainset_.size());
}

template<class Abc>
void CSClustApp<Abc>::InitLibrary() {
    if (!opts_.libfile.empty()) {  // read library from jumpstart file
        fprintf(out_, "Reading context library from %s ...\n",
                GetBasename(opts_.libfile).c_str());

        FILE* fin = fopen(opts_.libfile.c_str(), "r");
        if (!fin) throw Exception("Unable to read file '%s'!", opts_.libfile.c_str());
        lib_.reset(new ContextLibrary<Abc>(fin));
        fclose(fin);

    } else if (opts_.gauss_init != 0.0) {  // init by sampling weights from gaussian
        fputs("Initializing libray by sampling profiles from gaussian ...\n", out_);
        GaussianLibraryInit<Abc> init(opts_.gauss_init, *sm_, opts_.seed);
        lib_.reset(new ContextLibrary<Abc>(opts_.nprofs,
                                           trainset_.front().counts.length(),
                                           init));

    } else {  // init by sampling from training set
        fputs("Initializing library by sampling from training set ...\n", out_);
        MatrixPseudocounts<Abc> pc(*sm_);
        ConstantAdmix admix(opts_.pc_init);
        SamplingLibraryInit<Abc> init(trainset_, pc, admix, opts_.seed);
        lib_.reset(new ContextLibrary<Abc>(opts_.nprofs,
                                           trainset_.front().counts.length(),
                                           init));
    }
}

template<class Abc>
int CSClustApp<Abc>::Run() {
    InitSubstitutionMatrix();
    ReadTrainingData();
    InitLibrary();

    // Run EM clustering
    fprintf(out_, "Clustering %zu profile windows of length %zu into %zu profiles ...\n\n", trainset_.size(), lib_->wlen(), lib_->size());
    EMClustering<Abc> em(trainset_, *lib_, *sm_, opts_.em);
    RunClustering(em, opts_.toll, opts_.max_iters, out_);
    TransformToLin(em.lib);
    em.lib.SortByEntropy();

    // Copy learned library and add pseudocounts for color-space learning
    MatrixPseudocounts<Abc> pc(*sm_);
    ContextLibrary<Abc> lib_pc(em.lib);
    for (size_t k = 0; k < lib_pc.size(); ++k) {
        CountProfile<Abc> cp(lib_pc[k].probs);
        lib_pc[k].probs = pc.AddTo(cp, ConstantAdmix(0.5));
    }

    // Learn color-space SOM of context profiles
    const size_t ncolors = lib_pc.size() == AS219::kSize ? 4 : 5;
    GaussianLibraryInit<Abc> init(0.1, *sm_, opts_.seed);
    ContextLibrary<Abc> som(ncolors * ncolors * ncolors, lib_pc.wlen(), init);
    CoEmission<Abc> co_emission(lib_pc.wlen(), *sm_,
                                opts_.em.weight_center,
                                opts_.em.weight_decay);
    LearnContextMap(lib_pc, som, co_emission, 20000, (ncolors - 1) / 2, 0.1);

    // Map each context profile to an RGB color
    AssignContextColors(lib_pc, som, co_emission);

    // Assign names either by color or by abstract alphabet letter
    if (lib_pc.size() != AS219::kSize) {
        AssignContextNames(lib_pc, som, co_emission);
    } else {
        for (size_t k = 0; k < AS219::kSize; ++k)
            lib_pc[k].name = strprintf("%zu", k % 10);
    }

    // Copy names and colors from 'lib_pc' to real library WITHOUT pseudocounts
    for (size_t k = 0; k < lib_pc.size(); ++k) {
        em.lib[k].name = lib_pc[k].name;
        em.lib[k].color = lib_pc[k].color;
    }

    // Write profile library to outfile
    FILE* fout = fopen(opts_.outfile.c_str(), "w");
    if (!fout) throw Exception("Can't write outfile '%s'!", opts_.outfile.c_str());
    em.lib.Write(fout);
    fclose(fout);
    fprintf(out_, "\nWrote context library to %s\n", opts_.outfile.c_str());

    // Generate PDF with profile histograms
    if (!opts_.pdffile.empty()) {
        ContextLibraryPdfWriter<Abc> lib_writer(em.lib);
        if (em.lib.size() == AS219::kSize) {
            lib_writer.ncols = 73;
            lib_writer.margin = 0.5;
        } else if (em.lib.size() > 500 && Abc::kSize == AA::kSize) {
            lib_writer.pmin = 0.01;
        }
        lib_writer.WriteToFile(opts_.pdffile);
        fprintf(out_, "Wrote profile histograms to %s\n",  opts_.pdffile.c_str());
    }

    return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
    string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
    if (alphabet == "dna" || alphabet == "DNA")
        return cs::CSClustApp<cs::Dna>().main(argc, argv, stdout, "csclust");
    else
        return cs::CSClustApp<cs::AA>().main(argc, argv, stdout, "csclust");
}
