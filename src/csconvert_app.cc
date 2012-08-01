/*
  Copyright 2012 Christof Angermueller

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
#include "crf-inl.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CSConvertAppOptions {
  CSConvertAppOptions() { Init(); }
  virtual ~CSConvertAppOptions() {}

  // Set csconvert default parameters
  void Init() { 
    weight_center = 1.6;
    weight_decay  = 0.85;
    neff          = 1.0;
  }

  // Validates the parameter settings and throws exception if needed.
  void Validate() {
    if (infile.empty()) throw Exception("No input file provided!");
    if (neff < 1.0) throw Exception("Neff must be greater or equal 1.0!");
  }

  string infile;        // filename of the input model.
  string outfile;       // filename of the output model.
  double weight_center; // weight of the central profile column.
  double weight_decay;  // exponential decay of column weights.
  double neff;          // Neff sequences used for training the CRF.
};  // CSConvertAppOptions


template<class Abc>
class CSConvertApp : public Application {
 private:
  // Runs the count_profile_neff application.
  virtual int Run();
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& ops);
  // Prints options summary to stream.
  virtual void PrintOptions() const;
  // Prints short application description.
  virtual void PrintBanner() const;
  // Prints usage banner to stream.
  virtual void PrintUsage() const;
  // Converts a context library into a CRF.
  void Lib2Crf();
  // Converts a CRF into a context library.
  void Crf2Lib();

  // Parameter wrapper
  CSConvertAppOptions opts_;
};  // class CSConvertApp



template<class Abc>
void CSConvertApp<Abc>::ParseOptions(GetOpt_pp& ops) {
  ops >> Option('i', "infile", opts_.infile, opts_.infile);
  ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
  ops >> Option(' ', "weight-center", opts_.weight_center, opts_.weight_center);
  ops >> Option(' ', "weight-decay", opts_.weight_decay, opts_.weight_decay);
  ops >> Option(' ', "neff", opts_.neff, opts_.neff);
  opts_.Validate();
}

template<class Abc>
void CSConvertApp<Abc>::PrintBanner() const {
  fputs("Converts a context library into a CRF or vice versa.\n", out_);
}

template<class Abc>
void CSConvertApp<Abc>::PrintUsage() const {
  fputs("Usage: csconvert -i <infile> [options]\n", out_);
}

template<class Abc>
void CSConvertApp<Abc>::PrintOptions() const {
  fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
          "Input context library (*.lib) or CRF (*.crf)");
  fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>",
          "Filename of the output model");
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --weight-center [0;inf[",
          "Weight of the central profile column", opts_.weight_center);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --weight-decay [0;inf[",
          "Parameter for exponential decay of window weights", opts_.weight_decay);
  fprintf(out_, "  %-30s %s (def=%.2f)\n", "    --neff [1;inf[",
          "Mean Neff in sequences used for training the CRF", opts_.neff);
}

template<class Abc>
void CSConvertApp<Abc>::Lib2Crf() {
  // Read context library
  fputs("Reading context library ...\n", out_);
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (fin == NULL) throw Exception("Can't read from '%s'!", opts_.infile.c_str());
  ContextLibrary<Abc> lib(fin);
  fclose(fin);
  // Convert to CRF
  fputs("Converting context library to CRF ...\n", out_);
  LibraryBasedCrfInit<Abc> init(lib, opts_.weight_center, opts_.weight_decay);
  Crf<Abc> crf(lib.size(), lib.wlen(), init);
  // Write CRF
  fputs("Writing CRF ...\n", out_);
  string outfile = opts_.outfile;
  if (outfile.empty()) 
    outfile = GetBasename(opts_.infile, false) + ".crf";
  FILE* fout = fopen(outfile.c_str(), "w");
  if (fout == NULL) throw Exception("Can't write to '%s'!", outfile.c_str());
  crf.Write(fout);
  fclose(fout);
  fputs("Done!\n", out_);
}

template<class Abc>
void CSConvertApp<Abc>::Crf2Lib() {
  // Read context library
  fputs("Reading CRF ...\n", out_);
  FILE* fin = fopen(opts_.infile.c_str(), "r");
  if (fin == NULL) throw Exception("Can't read from '%s'!", opts_.infile.c_str());
  Crf<Abc> crf(fin);
  fclose(fin);
  // Convert to context library
  fputs("Converting CRF to context library ...\n", out_);
  CrfBasedLibraryInit<Abc> init(crf, opts_.weight_center, opts_.weight_decay, opts_.neff);
  ContextLibrary<Abc> lib(crf.size(), crf.wlen(), init);
  // Write context library
  fputs("Writing context library ...\n", out_);
  string outfile = opts_.outfile;
  if (outfile.empty()) 
    outfile = GetBasename(opts_.infile, false) + ".lib";
  FILE* fout = fopen(outfile.c_str(), "w");
  if (fout == NULL) throw Exception("Can't write to '%s'!", outfile.c_str());
  lib.Write(fout);
  fclose(fout);
  fputs("Done!\n", out_);
}

template<class Abc>
int CSConvertApp<Abc>::Run() {
  if (GetFileExt(opts_.infile) == "lib") Lib2Crf();
  else if (GetFileExt(opts_.infile) == "crf") Crf2Lib();
  else throw Exception("Invalid file extension!");
  return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
  string alphabet(getenv("CS_ALPHABET") ? getenv("CS_ALPHABET") : "");
  if (alphabet == "dna" || alphabet == "DNA")
    return cs::CSConvertApp<cs::Dna>().main(argc, argv, stdout, "csconvert");
  else
    return cs::CSConvertApp<cs::AA>().main(argc, argv, stdout, "csconvert");
}
