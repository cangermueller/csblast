#ifndef CS_CSBUILD_H
#define CS_CSBUILD_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation for HMM training.

#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "alignment.h"
#include "counts_profile.h"
#include "exception.h"
#include "getopt_pp.h"
#include "shared_ptr.h"
#include "utils.h"

using namespace GetOpt;

namespace cs
{

struct Params
{
    Params()
            : format("auto"),
              matchcol_assignment(-1),
              global_weights(false),
              log_level(Log::reporting_level())
    { }

    virtual ~Params() { }

    // Checks if all parameters are valid for running the training.
    void check();
    // Parses command line arguments and sets parameters accordingly.
    void parse_options(GetOpt_pp& options);

    // The input alignment file with training data.
    std::string infile;
    // The output file for the trained HMM.
    std::string outfile;
    // File format of input alignment
    std::string format;
    // Match column assignment for FASTA alignments
    int matchcol_assignment;
    // Use global instead of position specific weights for profile construction.
    bool global_weights;
    // The reporting level for logging.
    int log_level;
};



void Params::parse_options(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile, infile);
    options >> Option('o', "outfile", outfile, outfile);
    options >> Option('f', "format", format, format);
    options >> Option('M', "matchcol", matchcol_assignment, matchcol_assignment);
    options >> OptionPresent(' ', "global-weights", global_weights);
    options >> Option(' ', "log-level", log_level, log_level);
    Log::reporting_level() = Log::from_integer(log_level);

    check();
    if (outfile.empty()) outfile = get_file_basename(infile, false) + "prf";
    if (format == "auto") format = get_file_ext(infile);
}

void Params::check()
{
    if (infile.empty()) throw Exception("No input file provided!");
}



std::ostream& usage(const Params& params, std::ostream& out = std::cout)
{
    out << "Build a sequence profile from an alignment or sequence.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: csbuild -i <infile> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-30s %s\n",             "-i, --infile <filename>", "Path to input file with alignment or sequence");
    out << strprintf("  %-30s %s\n",             "-o, --outfile <filename>", "Path for output file with serialized profile");
    out << strprintf("  %-30s %s (def=%s)\n",    "-f, --format <string>", "Input data format: seq, fas, a2m, or a3m",
                     params.format.c_str());
    out << strprintf("  %-30s %s\n",             "-M, --matchcol [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    out << strprintf("  %-30s %s\n",             "", "(def: make columns with residue in first sequence match columns)");
    out << strprintf("  %-30s %s\n",             "    --global-weights", "Use global instead of position-specific weights for profiles");
    out << strprintf("  %-30s %s (def=%i)\n",    "    --log-level <int>", "Maximal reporting level for logging",
                     params.log_level);

    return out;
}



template<class Alphabet_T>
void csbuild(const Params& params, std::ostream& out = std::cout)
{
    std::ifstream fin(params.infile.c_str());
    if (!fin) throw Exception("Unable to read from input file '%s'!", params.infile.c_str());
    std::ofstream fout(params.outfile.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write profiles to output file '%s'!", params.outfile.c_str());

    if (params.format == "seq") {
        Sequence<Alphabet_T> seq(fin);
        CountsProfile<Alphabet_T> profile(seq);
        profile.write(fout);
        out << strprintf("Wrote profile with %i columns to %s\n", profile.num_cols(), params.outfile.c_str());
        LOG(INFO) << strprintf("Wrote profile with %i columns to %s\n", profile.num_cols(), params.outfile.c_str());

    } else {
        typename Alignment<Alphabet_T>::Format f = alignment_format_from_string<Alphabet_T>(params.format);
        Alignment<Alphabet_T> ali(fin, f);
        if (f == Alignment<Alphabet_T>::FASTA) {
            if (params.matchcol_assignment < 0)
                ali.assign_match_columns_by_sequence();
            else
                ali.assign_match_columns_by_gap_rule(params.matchcol_assignment);
        }
        CountsProfile<Alphabet_T> profile(ali, !params.global_weights);
        profile.write(fout);
        out << strprintf("Wrote profile with %i columns to %s\n", profile.num_cols(), params.outfile.c_str());
        LOG(INFO) << strprintf("Wrote profile with %i columns to %s\n", profile.num_cols(), params.outfile.c_str());
    }

    fout.close();
    fin.close();
}

}  // cs

#endif
