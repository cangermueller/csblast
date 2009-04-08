#ifndef CS_CSCLUST_H
#define CS_CSCLUST_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// CLI for EM clustering.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "alignment-inl.h"
#include "blosum_matrix.h"
#include "clustering.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "getopt_pp.h"
#include "matrix_pseudocounts.h"
#include "nucleotide_matrix.h"
#include "profile_library.h"
#include "shared_ptr.h"
#include "utils.h"

using namespace GetOpt;

namespace cs
{

struct Params : public ClusteringParams
{
    Params()
            : num_profiles(0),
              lib_pseudocounts(1.0f),
              data_pseudocounts(0.01f),
              blosum_type("BLOSUM62"),
              nucleotide_match(1.0f),
              nucleotide_mismatch(-3.0f),
              log_level(Log::reporting_level())
    { }

    virtual ~Params() { }

    // Checks if all parameters are valid for running the clustering.
    void check();
    // Parses command line arguments and sets parameters accordingly.
    void parse_options(GetOpt_pp& options);

#ifdef _WIN32
    static const char DIR_SEP = '\\';
#else
    static const char DIR_SEP = '/';
#endif

    // The input alignment file with training data.
    std::string infile;
    // The output file for profile library.
    std::string outfile;
    // Directory for output and temporary files
    std::string directory;
    // Library input file for restarting
    std::string libfile;
    // The number of profiles in the profile library.
    int num_profiles;
    // Pseudocounts to be added to each library profile.
    float lib_pseudocounts;
    // Pseudocounts to be added to observed data counts.
    float data_pseudocounts;
    // BLOSUM matrix for pseudocount generation.
    std::string blosum_type;
    // Reward for a nucleotide match.
    float nucleotide_match;
    // Penalty for a nucleotide mismatch
    float nucleotide_mismatch;
    // The reporting level for logging.
    int log_level;
};



void Params::parse_options(GetOpt_pp& options)
{
    options >> Option('i', "infile", infile, infile);
    options >> Option('o', "outfile", outfile, outfile);
    options >> Option('d', "directory", directory, directory);
    options >> Option('K', "num-profile", num_profiles, num_profiles);
    options >> Option('l', "likelihod-change", log_likelihood_change, log_likelihood_change);
    options >> Option('j', "jumpstart", libfile, libfile);
    options >> Option('B', "blocks", num_blocks, num_blocks);
    options >> Option('m', "matrix", blosum_type, blosum_type);
    options >> Option('q', "mismatch-score", nucleotide_mismatch, nucleotide_mismatch);
    options >> Option('r', "match-score", nucleotide_match, nucleotide_match);
    options >> Option(' ', "data-pc", data_pseudocounts, data_pseudocounts);
    options >> Option(' ', "lib-pc", lib_pseudocounts, lib_pseudocounts);
    options >> Option(' ', "min-scans", min_scans, min_scans);
    options >> Option(' ', "max-scans", max_scans, max_scans);
    options >> Option(' ', "weight-center", weight_center, weight_center);
    options >> Option(' ', "weight-decay", weight_decay, weight_decay);
    options >> Option(' ', "epsilon", epsilon_null, epsilon_null);
    options >> Option(' ', "beta", beta, beta);
    options >> Option(' ', "log-level", log_level, log_level);
    Log::reporting_level() = Log::from_integer(log_level);

    check();
    if (!directory.empty() && *directory.rbegin() != DIR_SEP) directory.append(1, DIR_SEP);
    if (outfile.empty()) outfile = directory + get_file_basename(infile, false) + "lib";
}

void Params::check()
{
    if (num_profiles == 0 && libfile.empty()) throw Exception("No value for number of profiles provided!");
    if (infile.empty()) throw Exception("No input file with training data provided!");
}



template<class Alphabet_T>
std::ostream& usage(const Params& params, std::ostream& out = std::cout);

template<class Alphabet_T>
std::ostream& substitution_matrix_options(const Params& params, std::ostream& out);

template<>
std::ostream& substitution_matrix_options<AminoAcid>(const Params& params, std::ostream& out);

template<class Alphabet_T>
std::ostream& usage(const Params& params, std::ostream& out = std::cout)
{
    out << "Cluster a training set of profile windows into a profile library.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: csclust -i <infile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-30s %s\n",             "-i, --infile <filename>", "Path to input file with profile windows");
    out << strprintf("  %-30s %s\n",             "-o, --outfile <filename>", "Path for output file with profile library");
    out << strprintf("  %-30s %s (def=%s)\n",    "-d, --directory <directory>", "Directory for temporary and output files",
                     params.directory.empty() ? "." : params.directory.c_str());
    out << strprintf("  %-30s %s\n",             "-K, --num-profiles [0,inf[", "Number of profiles in the library");
    out << strprintf("  %-30s %s (def=%3.1g)\n", "-l, --likelihood [0,inf[", "Maximal likelihood change per column for convergence",
                     params.log_likelihood_change);
    out << strprintf("  %-30s %s\n",             "-j, --jumpstart <filename>", "Jumpstart the clustering with a profile library.");
    out << strprintf("  %-30s %s\n",             "-B, --blocks [0,N]", "Number of blocks for online training (def: B=N^3/8)");

    substitution_matrix_options<Alphabet_T>(params, out);

    out << strprintf("  %-30s %s (def=%i)\n",    "    --min-scans [0,inf[", "Minimal number of training data scans",
                     params.min_scans);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --max-scans [0,inf[", "Maximal number of training data scans",
                     params.max_scans);
    out << strprintf("  %-30s %s (def=%3.1f)\n", "    --lib-pc [0,1]", "Pseudocounts for library profiles",
                     params.lib_pseudocounts);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --data-pc [0,1]", "Pseudocounts for training data",
                     params.data_pseudocounts);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --weight-center [0,1]", "Weight of central profile column in context window",
                     params.weight_center);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --weight-decay [0,1]", "Exponential decay of positional window weights",
                     params.weight_decay);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --epsilon [0,1]", "Start value for learning rate epsilon in online training",
                     params.epsilon_null);
    out << strprintf("  %-30s %s (def=%4.2f)\n", "    --beta [0,1]", "Exponential decay of epsilon in online training",
                     params.beta);
    out << strprintf("  %-30s %s (def=%i)\n",    "    --log-level <int>", "Maximal reporting level for logging",
                     params.log_level);

    return out;
}

template<class Alphabet_T>
std::ostream& substitution_matrix_options(const Params& params, std::ostream& out)
{
    out << strprintf("  %-30s %s (def=%.0f)\n",  "-q, --mismatch-score <int>", "Penalty for a nucleotide mismatch",
                     params.nucleotide_mismatch);
    out << strprintf("  %-30s %s (def=%.0f)\n",  "-r, --match-score <int>", "Reward for a nucleotide match",
                     params.nucleotide_match);
    return out;
}

template<>
std::ostream& substitution_matrix_options<AminoAcid>(const Params& params, std::ostream& out)
{
    out << strprintf("  %-30s %s (def=%s)\n",    "-m, --matrix <string>", "Substitution matrix: BLOSUM45, BLOSUM62, or BLOSUM80",
                     params.blosum_type.c_str());
    return out;
}



template<class Alphabet_T>
void csclust(const Params& params, std::ostream& out = std::cout);

template<class Alphabet_T>
shared_ptr< SubstitutionMatrix<Alphabet_T> > get_substitution_matrix(const Params& params);

template<>
shared_ptr< SubstitutionMatrix<AminoAcid> > get_substitution_matrix<AminoAcid>(const Params& params);

template<class Alphabet_T>
void csclust(const Params& params, std::ostream& out)
{
    typedef std::vector< shared_ptr< CountProfile<Alphabet_T> > > counts_vector;
    typedef typename counts_vector::iterator counts_iterator;

    shared_ptr< ProfileLibrary<Alphabet_T> > lib_ptr;
    shared_ptr< SubstitutionMatrix<Alphabet_T> > sm_ptr(get_substitution_matrix<Alphabet_T>(params));
    MatrixPseudocounts<Alphabet_T> sm_pc(sm_ptr.get());
    counts_vector data;

    // read data counts directly from serialized counts profiles
    std::ifstream fin(params.infile.c_str());
    if (!fin) throw Exception("Unable to read from input file '%s'!", params.infile.c_str());
    out << strprintf("Reading training profiles from %s ...", get_file_basename(params.infile).c_str());
    out.flush();
    LOG(INFO) << strprintf("Reading training profiles from %s ...", get_file_basename(params.infile).c_str());
    CountProfile<Alphabet_T>::readall(fin, data);
    out << strprintf(" %i profiles read", data.size()) << std::endl;
    LOG(INFO) << strprintf("%i profiles read", data.size());
    fin.close();

    // construct or read the profile library
    if (params.libfile.empty()) {
        out << strprintf("Initializing profile library by sampling %i profiles from training profiles ...", params.num_profiles);
        out.flush();
        LOG(INFO) << strprintf("Initializing HMM by sampling %i profiles from training profiles ...", params.num_profiles);
        SamplingProfileInitializer<Alphabet_T> profile_init(data, &sm_pc, params.lib_pseudocounts);
        lib_ptr = shared_ptr< ProfileLibrary<Alphabet_T> >(new ProfileLibrary<Alphabet_T>(params.num_profiles,
                                                                                          data[0]->num_cols(),
                                                                                          profile_init));
        lib_ptr->transform_to_logspace();
        out << std::endl;
    } else {
        std::ifstream fin(params.libfile.c_str());
        out << strprintf("Reading profile library from %s ...", get_file_basename(params.libfile).c_str());
        out.flush();
        LOG(INFO) << strprintf("Reading profile library from %s ...", get_file_basename(params.libfile).c_str());
        lib_ptr = shared_ptr< ProfileLibrary<Alphabet_T> >(new ProfileLibrary<Alphabet_T>(fin));
        out << std::endl;
    }

    // add pseudocounts to training data
    out << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", params.data_pseudocounts);
    out.flush();
    LOG(INFO) << strprintf("Adding pseudocounts to training profiles (admixture=%.2f) ...", params.data_pseudocounts);
    for (counts_iterator ci = data.begin(); ci != data.end(); ++ci) {
        sm_pc.add_to_profile(ConstantAdmixture(params.data_pseudocounts), ci->get());
        (*ci)->convert_to_counts();
    }
    out << std::endl;

    // cluster training data
    out << strprintf("Clustering training data by expectation maximization (K=%i, W=%i, N=%i) ...",
                     lib_ptr->num_profiles(), lib_ptr->num_cols(), data.size());
    out.flush();
    LOG(INFO) << strprintf("Clustering training data by expectation maximization (K=%i, W=%i, N=%i) ...",
                           lib_ptr->num_profiles(), lib_ptr->num_cols(), data.size());
    out << std::endl << std::endl;
    Clustering<Alphabet_T, CountProfile> em_clust(params, data, *lib_ptr, out);
    em_clust.run();

    // write profile library to outfile
    std::ofstream fout(params.outfile.c_str(), std::ios_base::out);
    if (!fout) throw Exception("Unable to write profile library to output file '%s'!", params.outfile.c_str());
    lib_ptr->write(fout);
    fout.close();
    out << std::endl << strprintf("Wrote profile library to %s", params.outfile.c_str()) << std::endl;
    LOG(INFO) << strprintf("Wrote profile_library to %s", params.outfile.c_str());
}

template<class Alphabet_T>
shared_ptr< SubstitutionMatrix<Alphabet_T> > get_substitution_matrix(const Params& params)
{
    return shared_ptr< SubstitutionMatrix<Alphabet_T> >(new NucleotideMatrix(params.nucleotide_match, params.nucleotide_mismatch));
}

template<>
shared_ptr< SubstitutionMatrix<AminoAcid> > get_substitution_matrix<AminoAcid>(const Params& params)
{
    BlosumMatrix::Type type = blosum_matrix_type_from_string(params.blosum_type);
    shared_ptr< SubstitutionMatrix<AminoAcid> > matrix_ptr(new BlosumMatrix(type));
    return matrix_ptr;
}


}  // cs

#endif
