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
#include "amino_acid.h"
#include "counts_profile.h"
#include "exception.h"
#include "getopt_pp.h"
#include "nucleotide.h"
#include "shared_ptr.h"
#include "util.h"

using namespace GetOpt;

namespace cs
{

struct Params
{
    Params()
            : format("auto"),
              num_profiles(std::numeric_limits<int>::max()),
              random_shuffle(false),
              window_length(0),
              sample_rate(0.2f),
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
    // The number of profiles to be converted.
    int num_profiles;
    // Randomly reorder the input alignments/sequences?
    bool random_shuffle;
    // The number of columns in the context window.
    int window_length;
    // Fraction of profile windows sampled from each full length alignment/sequence.
    float sample_rate;
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
    options >> Option('N', "num-profiles", num_profiles, num_profiles);
    options >> OptionPresent('r', "random-shuffle", random_shuffle);
    options >> Option('W', "window-length", window_length, window_length);
    options >> Option('s', "sample-rate", sample_rate, sample_rate);
    options >> Option('M', "matchcol-assignment", matchcol_assignment, matchcol_assignment);
    options >> OptionPresent(' ', "global-weights", global_weights);
    options >> Option(' ', "log-level", log_level, log_level);
    Log::reporting_level() = Log::from_integer(log_level);

    check();
    if (outfile.empty()) outfile = get_file_basename(infile, false) + "prf";
    if (format == "auto") format = get_file_ext(infile);
}

void Params::check()
{
    if (infile.empty()) throw Exception("No input file with training data provided!");
}



std::ostream& usage(const Params& params, std::ostream& out = std::cout)
{
    out << "Build a profile (database) with counts calculated from an alignment or sequence (database).\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: csbuild -i <infile> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-38s %s\n",             "-i, --infile <filename>", "Path to input file with alignments or sequences");
    out << strprintf("  %-38s %s\n",             "-o, --outfile <filename>", "Path to output file with profiles");
    out << strprintf("  %-38s %s (def=%s)\n",    "-f, --format <string>", "Input data format: seq, fas, a2m, or a3m",
                     params.format.c_str());
    out << strprintf("  %-38s %s\n",             "-N, --num-profiles [0,inf[", "Limit profile construction to N profiles");
    out << strprintf("  %-38s %s\n",             "-r, --random-shuffle", "Randomly reorder input alignments/sequences before conversion.");
    out << strprintf("  %-38s %s\n",             "-W, --window-length [0,inf[", "Extract profile windows of length W instead of full-length profiles");
    out << strprintf("  %-38s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]", "Fraction of profile windows sampled from each alignment/sequence",
                     params.sample_rate);
    out << strprintf("  %-38s %s\n",             "-M, --matchcol-assignment [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    out << strprintf("  %-38s %s\n",             "", "(def: make columns with residue in first sequence match columns)");
    out << strprintf("  %-38s %s\n",             "    --global-weights", "Use global instead of position-specific weights for profiles");
    out << strprintf("  %-38s %s (def=%i)\n",    "    --log-level <int>", "Maximal reporting level for logging",
                     params.log_level);

    return out;
}



template<class Alphabet_T>
void csbuild(const Params& params, std::ostream& out = std::cout);

template< class Alphabet_T,
          template<class Alphabet_U> class Input_T >
void build_profiles(const Params& params,
                    const std::vector< shared_ptr< Input_T<Alphabet_T> > >& data,
                    std::vector< shared_ptr< CountsProfile<Alphabet_T> > >& profiles,
                    std::ostream& out);

template<class Alphabet_T>
shared_ptr< CountsProfile<Alphabet_T> > create_profile(const Alignment<Alphabet_T>& ali, const Params& params);

template<class Alphabet_T>
shared_ptr< CountsProfile<Alphabet_T> > create_profile(const Sequence<Alphabet_T>& seq, const Params& params);

template<class Alphabet_T>
void csbuild(const Params& params, std::ostream& out)
{
    typedef std::vector< shared_ptr< CountsProfile<Alphabet_T> > > profiles_vector;
    typedef typename profiles_vector::const_iterator profile_iterator;

    std::vector< shared_ptr< CountsProfile<Alphabet_T> > > profiles;
    std::ifstream fin(params.infile.c_str());

    // read input data from infile
    if (params.format == "seq") {
        // read sequences
        out << strprintf("Reading sequences from %s ...\n", get_file_basename(params.infile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading sequences from %s ...", get_file_basename(params.infile).c_str());

        std::vector< shared_ptr<Sequence<Alphabet_T> > > seqs;
        int i = 0;
        while (fin.peek() && fin.good()) {
            shared_ptr< Sequence<Alphabet_T> > seq_ptr(new Sequence<Alphabet_T>(fin));
            seqs.push_back(seq_ptr);

            i += 1;
            if (i % 2 == 0) {
                out << '.';
                out.flush();
            }
            if (i % 100 == 0) out << " " << i << std::endl;
        }
        if (i % 100 != 0) out << std::string(50 - iround((i % 100) / 2), ' ') << " " << i << std::endl;

        // build profiles
        if (params.random_shuffle) random_shuffle(seqs.begin(), seqs.end());
        build_profiles(params, seqs, profiles, out);

    } else {
        // read alignments
        out << strprintf("Reading alignments from %s ...\n", get_file_basename(params.infile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading sequences from %s ...", get_file_basename(params.infile).c_str());

        typename Alignment<Alphabet_T>::Format f = alignment_format_from_string<Alphabet_T>(get_file_ext(params.infile));
        std::vector< shared_ptr<Alignment<Alphabet_T> > > alis;
        int i = 0;
        while (fin.peek() && fin.good()) {
            shared_ptr< Alignment<Alphabet_T> > ali_ptr(new Alignment<Alphabet_T>(fin, f));
            if (f == Alignment<Alphabet_T>::FASTA) {
                if (params.matchcol_assignment < 0)
                    ali_ptr->assign_match_columns_by_sequence();
                else
                    ali_ptr->assign_match_columns_by_gap_rule(params.matchcol_assignment);
            }
            ali_ptr->remove_insert_columns();  // remove insert columns to save memory!
            alis.push_back(ali_ptr);

            i += 1;
            if (i % 2 == 0) {
                out << '.';
                out.flush();
            }
            if (i % 100 == 0) out << " " << i << std::endl;
        }
        if (i % 100 != 0) out << std::string(50 - iround((i % 100) / 2), ' ') << " " << i << std::endl;

        // build profiles
        if (params.random_shuffle) random_shuffle(alis.begin(), alis.end());
        build_profiles(params, alis, profiles, out);
    }

    // write HMM to outfile
    std::fstream fout(params.outfile.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write profiles to output file '%s'!", params.outfile.c_str());

    int num_cols = 0;
    for (profile_iterator it = profiles.begin(); it != profiles.end(); ++it) {
        (*it)->write(fout);
        num_cols += (*it)->num_cols();
    }

    out << std::endl << strprintf("Wrote %i profiles with a total number of %i columns to %s\n",
                                  profiles.size(), num_cols, params.outfile.c_str());
    LOG(INFO) << strprintf("Wrote %i profiles with a total number of %i columns to %s",
                                  profiles.size(), num_cols, params.outfile.c_str());
    fout.close();
    fin.close();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Input_T >
void build_profiles(const Params& params,
                    const std::vector< shared_ptr< Input_T<Alphabet_T> > >& data,
                    std::vector< shared_ptr< CountsProfile<Alphabet_T> > >& profiles,
                    std::ostream& out)
{
    typedef typename std::vector< shared_ptr< Input_T<Alphabet_T> > >::const_iterator data_iterator;
    typedef typename std::vector<int>::const_iterator index_iterator;

    out << "Building profiles ...\n";
    out.flush();
    LOG(INFO) << "Building profiles ...";

    // iterate over input data and build counts profiles either by full-length conversion or by sampling of context windows
    int i = 0;
    for (data_iterator it = data.begin(); it != data.end() && static_cast<int>(profiles.size()) < params.num_profiles; ++it) {
        shared_ptr< CountsProfile<Alphabet_T> > profile_ptr = create_profile(**it, params);

        if (params.window_length == 0) {  // full length conversion
            profiles.push_back(profile_ptr);

            i += 1;
            if (i % 2 == 0) {
                out << '.';
                out.flush();
            }
            if (i % 100 == 0) out << " " << i << std::endl;

        } else if (profile_ptr->num_cols() >= params.window_length) {  // sampling of context windows
            // prepare sample of indices
            std::vector<int> idx;
            for (int j = 0; j <= profile_ptr->num_cols() - params.window_length; ++j) idx.push_back(j);
            random_shuffle(idx.begin(), idx.end());
            const int sample_size = iround(params.sample_rate * idx.size());
            idx.erase(idx.begin() + sample_size, idx.end());  // sample only a fraction of the profile indices.

            // Add sub-profiles at sampled indices to HMM
            for (index_iterator ii = idx.begin(); ii != idx.end() && static_cast<int>(profiles.size()) < params.num_profiles; ++ii) {
                shared_ptr< CountsProfile<Alphabet_T> > p(new CountsProfile<Alphabet_T>(*profile_ptr, *ii, params.window_length));
                profiles.push_back(p);

                i += 1;
                if (i % 2 == 0) {
                    out << '.';
                    out.flush();
                }
                if (i % 100 == 0) out << " " << i << std::endl;
            }
        }
    }
    if (i % 100 != 0) out << std::string(50 - iround((i % 100) / 2), ' ') << " " << i << std::endl;
}

template<class Alphabet_T>
inline shared_ptr< CountsProfile<Alphabet_T> > create_profile(const Alignment<Alphabet_T>& ali, const Params& params)
{
    return shared_ptr< CountsProfile<Alphabet_T> >(new CountsProfile<Alphabet_T>(ali, !params.global_weights));
}

template<class Alphabet_T>
inline shared_ptr< CountsProfile<Alphabet_T> > create_profile(const Sequence<Alphabet_T>& seq, const Params&)
{
    return shared_ptr< CountsProfile<Alphabet_T> >(new CountsProfile<Alphabet_T>(seq));
}

}  // cs

#endif
