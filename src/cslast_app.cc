/*
  Copyright 2011 Andreas Biegert

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
#include "context_library.h"
#include "crf_pseudocounts-inl.h"
#include "crf-inl.h"
#include "cslast.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "last_pssm.h"
#include "sequence-inl.h"
#include "tamura_nei_matrix.h"

using namespace GetOpt;
using std::string;
using std::map;
using std::vector;

namespace cs {

typedef vector<Sequence<Dna> > SeqVec;

struct CSLastAppOptions {

    CSLastAppOptions() { Init(); }

    void Init() {
        pc_admix        = 0.90;
        pc_ali          = 12.0;
        pc_engine       = "auto";
        weight_center   = 1.6;
        weight_decay    = 0.85;
    }

    // Validates parameter settings and throws exception if needed.
    void Validate() {
        if (infile.empty()) throw Exception("No input file provided!");
        if (modelfile.empty()) throw Exception("No context data provided!");
    }

    // The input query sequence.
    string infile;
    // The output file for results.
    string outfile;
    // Input file with context profile library or HMM
    string modelfile;
    // Database file
    string dbfile;
    // Overall pseudocount admixture
    double pc_admix;
    // Constant in pseudocount calculation for alignments
    double pc_ali;
    // Pseudocount engine
    string pc_engine;
    // Path to LAST executable
    string last_path;
    // Weight of central column in multinomial emission
    double weight_center;
    // Exponential decay of window weights
    double weight_decay;
    // LAST options map
    CSLastOptions cslast;
};  // struct CSLastAppOptions


class CSLastApp : public Application {
  private:
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
    // Initializes all class members for CSI-LAST searches
    void Init();
    // Writes current PSSM in LAST format
    void SavePssm() const;
    // Add pseudocounts and prepares CS-BLAST engine for run with given query
    void PrepareForRun(const Sequence<Dna>& query);

    // Parameter wrapper
    CSLastAppOptions opts_;
    // Query sequence
    scoped_ptr<Sequence<Dna> > query_;
    // Profile library for pseudocounts
    scoped_ptr<ContextLibrary<Dna> > lib_;
    // CRF for pseudocounts
    scoped_ptr<Crf<Dna> > crf_;
    // Pseudocount engine
    scoped_ptr<Pseudocounts<Dna> > pc_;
    // CS-BLAST engine
    scoped_ptr<CSLast> cslast_;
    // PSSM for LAST jumpstarting
    scoped_ptr<LastPssm> pssm_;
    // Substitution matrix for background frequencies
    scoped_ptr<SubstitutionMatrix<Dna> > mat_;
    // Vector with pointers to query sequences
    SeqVec queries_;
};  // class CSLastApp



void CSLastApp::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
    ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
    ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
    ops >> Option(' ', "pc-engine", opts_.pc_engine, opts_.pc_engine);
    ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
    ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
    ops >> Option(' ', "last-path", opts_.last_path, opts_.last_path);
    ops >> Option(' ', "LAST_PATH", opts_.last_path, opts_.last_path);

    vector<string> dbfile_infile;
    ops >> Option(GetOpt_pp::EMPTY_OPTION, dbfile_infile);
    if (dbfile_infile.size() > 2) {
        opts_.dbfile = dbfile_infile[0];
        opts_.infile = dbfile_infile[1];
        opts_.modelfile = dbfile_infile[2];
    }

    // Put remaining arguments into PSI-LAST options map
    for(GetOpt_pp::short_iterator it = ops.begin(); it != ops.end(); ++it) {
        if (!it.extracted())
            opts_.cslast[it.option()] = it.args().front();
    }

    opts_.Validate();

    if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
        opts_.pc_engine = GetFileExt(opts_.modelfile);
}

void CSLastApp::PrintBanner() const {
    fputs("Search with a DNA sequence against DNA databases for locally\n"
          "similar sequences.\n", out_);
}

void CSLastApp::PrintUsage() const {
    fputs("Usage: cslast <db-name> <sequence-file> <context-file> --last-path <lastdir> [options]\n", out_);
}

void CSLastApp::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
            "Output file with search results (def=stdout)");
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
            "Pseudocount admix for context-specific pseudocounts", opts_.pc_admix);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-center [0,inf[",
            "Weight of central profile column", opts_.weight_center);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "    --weight-decay [0,inf[",
            "Parameter for exponential decay of window weights", opts_.weight_decay);
    fprintf(out_, "  %-30s %s\n", "    --last-path <path>",
            "Path to directory with lastal executable (or set LAST_PATH)");
}

int CSLastApp::Run() {
    int status = 0;
    Init();

    for (SeqVec::iterator it = queries_.begin(); it != queries_.end(); ++it) {
        PrepareForRun(*it);

        cslast_->set_options(opts_.cslast);

        // Run CS-LAST
        FILE* fout = opts_.outfile.empty() ? out_ : fopen(opts_.outfile.c_str(), it != queries_.begin() ? "a" : "w");
        if (!fout) throw Exception("Unable to write to '%s'!", opts_.outfile.c_str());
        status = cslast_->Run(fout);
        if (!opts_.outfile.empty()) fclose(fout);
    }

    return status;
}

void CSLastApp::Init() {
    // Read query sequences
    FILE* fin = fopen(opts_.infile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!", opts_.infile.c_str());
    ReadAll(fin, queries_);
    fclose(fin);
    if (queries_.empty())
        throw Exception("No sequences found in '%s'!", opts_.infile.c_str());

    // Setup pseudocount engine
    if (opts_.pc_engine == "lib") {
        fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        lib_.reset(new ContextLibrary<Dna>(fin));
        fclose(fin);

        TransformToLog(*lib_);
        pc_.reset(new LibraryPseudocounts<Dna>(*lib_,
                                              opts_.weight_center,
                                              opts_.weight_decay));

    } else if (opts_.pc_engine == "crf") {
        fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin) throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        crf_.reset(new Crf<Dna>(fin));
        fclose(fin);

        pc_.reset(new CrfPseudocounts<Dna>(*crf_));

    } else {
        throw Exception("Unknown pseudocount engine '%s'!", opts_.pc_engine.c_str());
    }

    mat_.reset(new TamuraNeiMatrix());
}

void CSLastApp::PrepareForRun(const Sequence<Dna>& query) {
    // Setup PSSM of query profile with context-specific pseudocounts
    ConstantAdmix admix(opts_.pc_admix);
    pssm_.reset(new LastPssm(query, pc_->AddTo(query, admix), *mat_));

    // Setup CS-LAST engine
    cslast_.reset(new CSLast(opts_.dbfile, pssm_.get(), opts_.cslast));

    // Set path to PSI-BLAST executable
    if (!opts_.last_path.empty())
        cslast_->set_exec_path(opts_.last_path);
}

}  // namespace cs

int main(int argc, char* argv[]) {
    return cs::CSLastApp().main(argc, argv, stdout, "cslast");
}
