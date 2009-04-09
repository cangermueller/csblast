#ifndef CS_CSSAMPLE_H
#define CS_CSSAMPLE_H
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

#include "alignment-inl.h"
#include "count_profile-inl.h"
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
            : sample_size(std::numeric_limits<int>::max()),
              window_length(0),
              sample_rate(0.2f),
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
    // The number of profiles to be sampled out.
    int sample_size;
    // The number of columns in the context window.
    int window_length;
    // Fraction of profile windows sampled from each full length alignment/sequence.
    float sample_rate;
    // The reporting level for logging.
    int log_level;
};



void Params::parse_options(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile, infile);
    options >> Option('o', "outfile", outfile, outfile);
    options >> Option('N', "sample-size", sample_size, sample_size);
    options >> Option('W', "window-length", window_length, window_length);
    options >> Option('s', "sample-rate", sample_rate, sample_rate);
    options >> Option(' ', "log-level", log_level, log_level);
    Log::reporting_level() = Log::from_integer(log_level);

    check();
    if (outfile.empty()) outfile = get_file_basename(infile, false) + "prf";
}

void Params::check()
{
    if (infile.empty()) throw Exception("No input file provided!");
    if (outfile.empty()) throw Exception("No output file provided!");
}



std::ostream& usage(const Params& params, std::ostream& out = std::cout)
{
    out << "Sample (context) profiles from a large profile database.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cssample -i <infile> -o <outfile> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-30s %s\n",             "-i, --infile <filename>", "Path to input file with profile database");
    out << strprintf("  %-30s %s\n",             "-o, --outfile <filename>", "Path for output file with sampled profiles");
    out << strprintf("  %-30s %s\n",             "-N, --num-profiles [0,inf[", "Maximal number of profiles to sample (def=inf)");
    out << strprintf("  %-30s %s\n",             "-W, --window-length [0,inf[", "Sample context profiles of length W instead of full-length profiles");
    out << strprintf("  %-30s %s (def=%3.1f)\n", "-s, --sample-rate [0,1]", "Fraction of context profiles sampled per full-length profile",
                     params.sample_rate);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --log-level <int>", "Maximal reporting level for logging",
                     params.log_level);

    return out;
}



template<class Alphabet_T>
void cssample(const Params& params, std::ostream& out = std::cout);

template< class Alphabet_T,
          template<class Alphabet_U> class Input_T >
void sample_profiles(const Params& params,
                    const std::vector< shared_ptr< CountProfile<Alphabet_T> > >& profiles,
                    std::vector< shared_ptr< CountProfile<Alphabet_T> > >& samples,
                    std::ostream& out);

template<class Alphabet_T>
void cssample(const Params& params, std::ostream& out)
{
    typedef std::vector< shared_ptr< CountProfile<Alphabet_T> > > profiles_vector;
    typedef typename profiles_vector::const_iterator profile_iterator;

    std::vector< shared_ptr< CountProfile<Alphabet_T> > > profiles;  // large pool of profiles to sample from
    std::vector< shared_ptr< CountProfile<Alphabet_T> > > samples;   // sampled profiles

    // read database of profiles
    std::ifstream fin(params.infile.c_str());
    if (!fin) throw Exception("Unable to read from input file '%s'!", params.infile.c_str());
    out << strprintf("Reading profiles from %s ...", get_file_basename(params.infile).c_str());
    out.flush();
    LOG(INFO) << strprintf("Reading profiles from %s ...", get_file_basename(params.infile).c_str());
    CountProfile<Alphabet_T>::readall(fin, profiles);
    out << strprintf(" %i profiles read", profiles.size()) << std::endl;
    LOG(INFO) << strprintf("%i profiles read", profiles.size());
    fin.close();

    // sample profiles
    random_shuffle(profiles.begin(), profiles.end());
    sample_profiles(params, profiles, samples, out);

    // write sampled profiles to outfile
    std::ofstream fout(params.outfile.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write profiles to output file '%s'!", params.outfile.c_str());
    int num_cols = 0;
    for (profile_iterator it = samples.begin(); it != samples.end(); ++it) {
        (*it)->write(fout);
        num_cols += (*it)->num_cols();
    }
    out << strprintf("Wrote %i profiles with a total number of %i columns to %s\n",
                                  samples.size(), num_cols, params.outfile.c_str());
    LOG(INFO) << strprintf("Wrote %i profiles with a total number of %i columns to %s",
                                  samples.size(), num_cols, params.outfile.c_str());
    fout.close();
}

template<class Alphabet_T>
void sample_profiles(const Params& params,
                     const std::vector< shared_ptr< CountProfile<Alphabet_T> > >& profiles,
                     std::vector< shared_ptr< CountProfile<Alphabet_T> > >& samples,
                     std::ostream& out)
{
    typedef typename std::vector< shared_ptr< CountProfile<Alphabet_T> > >::const_iterator profile_iterator;
    typedef typename std::vector<int>::const_iterator index_iterator;

    out << strprintf("Sampling %i profiles from pool of %i profiles ...\n", params.sample_size, profiles.size());
    out.flush();
    LOG(INFO) << strprintf("Sampling %i profiles from pool of %i profiles ...\n", params.sample_size, profiles.size());

    // iterate over input data and build counts profiles either by full-length conversion or by sampling of context windows
    for (profile_iterator it = profiles.begin(); it != profiles.end() && static_cast<int>(samples.size()) < params.sample_size; ++it) {
        if (params.window_length == 0) {  // add full length profile to samples
            samples.push_back(*it);
        } else if ((*it)->num_cols() >= params.window_length) {  // sample context windows if profile has sufficient length
            // prepare sample of indices
            std::vector<int> idx;
            for (int j = 0; j <= (*it)->num_cols() - params.window_length; ++j) idx.push_back(j);
            random_shuffle(idx.begin(), idx.end());
            const int sample_size = iround(params.sample_rate * idx.size());
            idx.erase(idx.begin() + sample_size, idx.end());  // sample only a fraction of the profile indices.

            // Add sub-profiles at sampled indices to HMM
            for (index_iterator ii = idx.begin(); ii != idx.end() && static_cast<int>(samples.size()) < params.sample_size; ++ii) {
                shared_ptr< CountProfile<Alphabet_T> > p(new CountProfile<Alphabet_T>(**it, *ii, params.window_length));
                samples.push_back(p);
            }
        }
    }
}

}  // namespace cs

#endif
