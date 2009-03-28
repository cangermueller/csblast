#ifndef CS_CLUSTERING_H
#define CS_CLUSTERING_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation of expectation-maximization clustering.

#include <cmath>

#include <iostream>
#include <valarray>
#include <vector>

#include "context_profile.h"
#include "expectation_maximization.h"
#include "log.h"
#include "profile_matcher.h"
#include "profile_library.h"
#include "progress_table.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

// Forward declaration;
template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ClusteringProgressTable;


struct ClusteringParams : public ProfileMatcherParams, public ExpectationMaximizationParams
{
    ClusteringParams()
            : ProfileMatcherParams(),
              ExpectationMaximizationParams()
    { }

    ClusteringParams(const ClusteringParams& params)
            : ProfileMatcherParams(params),
              ExpectationMaximizationParams(params)
    { }
};


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class Clustering : public ExpectationMaximization<Alphabet_T, Subject_T>
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_vector;

    // Initializes a new clustering object without output.
    Clustering(const ClusteringParams& params,
               const data_vector& data,
               ProfileLibrary<Alphabet_T>& lib);
    // Initializes a new clustering object with output.
    Clustering(const ClusteringParams& params,
               const data_vector& data,
               ProfileLibrary<Alphabet_T>& lib,
               std::ostream& out);

    virtual ~Clustering();

  protected:
    // Needed to access names in templatized base class
    using ExpectationMaximization<Alphabet_T, Subject_T>::progress_table_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::log_likelihood_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::num_eff_cols_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::epsilon_;
    using ExpectationMaximization<Alphabet_T, Subject_T>::data_;

    // Evaluates the responsibilities using the current parameter values.
    virtual void expectation_step(const data_vector& block);
    // Reestimate teh parameters using the current responsibilities.
    virtual void maximization_step();
    // Prepares all members for clustering.
    virtual void init();
    // Returns parameter wrapper
    virtual const ClusteringParams& params() const { return params_; }
    // Adds the contribution of the responsibilities for a subject to sufficient statistics for priors.
    void add_contribution_to_priors(const std::valarray<double>& p_zn);
    // Adds the contribution of the responsibilities for a counts profile to sufficient statistics for emissions.
    void add_contribution_to_emissions(const std::valarray<double>& p_zn, const CountsProfile<Alphabet_T>& c);
    // Adds the contribution of the responsibilities for a sequence to sufficient statistics for emissions.
    void add_contribution_to_emissions(const std::valarray<double>& p_zn, const Sequence<Alphabet_T>& s);
    // Updates global sufficient statistics with sufficient statistics calculated on current block.
    void update_sufficient_statistics();

    // Parameter wrapper for clustering.
    const ClusteringParams& params_;
    // Profile library with context profiles
    ProfileLibrary<Alphabet_T> lib_;
    // Profile matcher for calculation of emission probabilities.
    ProfileMatcher<Alphabet_T> profile_matcher_;
    // Global expected sufficient statistics for emissions and state priors.
    profiles_vector profile_stats_;
    // Expeted sufficient statistics for emissions and state priors based on current block.
    profiles_vector profile_stats_block_;
};


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ClusteringProgressTable : public ProgressTable
{
  public:
    ClusteringProgressTable(const Clustering<Alphabet_T, Subject_T>* clustering,
                            std::ostream& out = std::cout,
                            int width = 30);

    virtual ~ClusteringProgressTable() { }

    virtual void print_header();
    virtual void print_row_begin();
    virtual void print_row_end();

  protected:
    const Clustering<Alphabet_T, Subject_T>* clustering_;
};



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
Clustering<Alphabet_T, Subject_T>::Clustering(const ClusteringParams& params,
                                                     const data_vector& data,
                                                     ProfileLibrary<Alphabet_T>& lib)
        : ExpectationMaximization<Alphabet_T, Subject_T>(data),
          params_(params),
          lib_(lib),
          profile_matcher_(lib.num_cols(), params),
          profile_stats_(),
          profile_stats_block_()
{
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
Clustering<Alphabet_T, Subject_T>::Clustering(const ClusteringParams& params,
                                                     const data_vector& data,
                                                     ProfileLibrary<Alphabet_T>& lib,
                                                     std::ostream& out)
        : ExpectationMaximization<Alphabet_T, Subject_T>(data),
          params_(params),
          lib_(lib),
          profile_matcher_(lib.num_cols(), params),
          profile_stats_(),
          profile_stats_block_()
{
    progress_table_ = new ClusteringProgressTable<Alphabet_T, Subject_T>(this, out);
    init();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
Clustering<Alphabet_T, Subject_T>::~Clustering()
{
    if (progress_table_) delete progress_table_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::expectation_step(const data_vector& block)
{
    const int num_profiles = lib_.num_profiles();
    std::valarray<double> p_zn(0.0f, lib_.num_profiles());

    // Given each training window compute posterior probabilities p_zn[k] of profile k
    for (typename data_vector::const_iterator bi = block.begin(); bi != block.end(); ++bi) {
        double sum = 0.0f;
        for (int k = 0; k < num_profiles; ++k) {
            p_zn[k] = lib_[k].prior() * profile_matcher_(lib_[k], **bi, lib_[k].center());
            sum += p_zn[k];

            LOG(DEBUG2) << strprintf("alpha(%i)=%-8.5g   P(c_n|p_%i)=%-8.5g   P(z_n=%-4i)=%-8.5g",
                                     k, lib_[k].prior(), k, profile_matcher_(lib_[k], **bi, lib_[k].center()), k, p_zn[k]);
        }
        p_zn /= sum;
        add_contribution_to_priors(p_zn);
        add_contribution_to_emissions(p_zn, **bi);
        log_likelihood_ += log(sum) / num_eff_cols_;
        LOG(DEBUG1) << strprintf("log(L)=%-8.5g", log_likelihood_);

        if (progress_table_) progress_table_->print_progress(lib_.num_profiles());
    }

    update_sufficient_statistics();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::maximization_step()
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    float sum = 0.0f;
    for (int k = 0; k < num_profiles; ++k) sum += profile_stats_[k]->prior();
    float fac = 1.0f / sum;

    // Assign new priors and emission probabilities
    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_[k];
        LOG(DEBUG1) << p_k;

        lib_[k].set_prior(p_k.prior() * fac);
        ContextProfile<Alphabet_T> tmp(p_k);
        if (normalize(tmp)) {  // don't update profiles that did'n get any evidence
            tmp.transform_to_logspace();
            for (int i = 0; i < num_cols; ++i)
                for (int a = 0; a < alphabet_size; ++a)
                    lib_[k][i][a] = tmp[i][a];
        }
    }

    ++lib_;  // Increment iteration counter
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_priors(const std::valarray<double>& p_zn)
{
    const int num_profiles = lib_.num_profiles();

    for (int k = 0; k < num_profiles; ++k)
        profile_stats_block_[k]->set_prior(profile_stats_block_[k]->prior() + p_zn[k]);
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_emissions(const std::valarray<double>& p_zn,
                                                                             const CountsProfile<Alphabet_T>& c)
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];

        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p_k[j][a] += c[j][a] * p_zn[k];
            }
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::add_contribution_to_emissions(const std::valarray<double>& p_zn,
                                                                             const Sequence<Alphabet_T>& s)
{
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_k = *profile_stats_block_[k];

        for (int j = 0; j < num_cols; ++j) {
            p_k[j][s[j]] += p_zn[k];
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline void Clustering<Alphabet_T, Subject_T>::update_sufficient_statistics()
{
    const float gamma       = 1.0f - epsilon_;
    const int num_profiles  = lib_.num_profiles();
    const int num_cols      = lib_.num_cols();
    const int alphabet_size = lib_.alphabet_size();

    // Update priors and emissions statistics
    for (int k = 0; k < num_profiles; ++k) {
        ContextProfile<Alphabet_T>& p_block = *profile_stats_block_[k];
        ContextProfile<Alphabet_T>& p       = *profile_stats_[k];

        p.set_prior(p.prior() * gamma + p_block.prior());
        for (int j = 0; j < num_cols; ++j) {
            for (int a = 0; a < alphabet_size; ++a) {
                p[j][a] = gamma * p[j][a] + p_block[j][a];
            }
        }
        reset(p_block);
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void Clustering<Alphabet_T, Subject_T>::init()
{
    // Create profiles for global and block-level sufficient statistics
    for (int k = 0; k < lib_.num_profiles(); ++k) {
        profile_stats_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, lib_.num_cols())));
        profile_stats_block_.push_back(shared_ptr< ContextProfile<Alphabet_T> >(new ContextProfile<Alphabet_T>(k, lib_.num_cols())));
    }

    // Initialize total amount of work per scan
    if (progress_table_) progress_table_->set_total_work(lib_.num_profiles() * data_.size());

    // Compute effective number of training data columns
    num_eff_cols_ = profile_matcher_.num_eff_cols() * data_.size();
}



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
ClusteringProgressTable<Alphabet_T, Subject_T>::ClusteringProgressTable(const Clustering<Alphabet_T, Subject_T>* clustering,
                                                                        std::ostream& out,
                                                                        int width)
        : ProgressTable(out, width),
          clustering_(clustering)
{ }

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ClusteringProgressTable<Alphabet_T, Subject_T>::print_header()
{
    out_ << strprintf("%-4s %4s %4s %7s  %-30s  %9s  %8s\n",
                      "Scan", "Itrs", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
    out_ << std::string(75, '-') << std::endl;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ClusteringProgressTable<Alphabet_T, Subject_T>::print_row_begin()
{
    reset();
    out_ << strprintf("%-4i %4i %4i %7.4f  ", clustering_->scan(), clustering_->iterations(),
                      clustering_->num_blocks(), clustering_->epsilon());
    out_.flush();
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ClusteringProgressTable<Alphabet_T, Subject_T>::print_row_end()
{
    if (clustering_->scan() == 1)
        out_ << strprintf("  %9.5f\n", clustering_->log_likelihood());
    else
        out_ << strprintf("  %9.5f  %+8.5f\n", clustering_->log_likelihood(), clustering_->log_likelihood_change());
    out_.flush();
}

}  // cs

#endif
