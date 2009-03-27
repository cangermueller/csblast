#ifndef CS_PROFILE_MATCHER_H
#define CS_PROFILE_MATCHER_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Encapsulation for matching a context profile with a counts profile or a
// sequence using multinomial distributions of alphabet probabilities.

#include <cmath>
#include <valarray>

#include "exception.h"
#include "context_profile.h"
#include "counts_profile.h"
#include "sequence.h"

namespace cs
{

struct ProfileMatcherParams
{
    ProfileMatcherParams()
            : ignore_context(false),
              weight_center(1.6f),
              weight_decay(0.85f)
    { }

    ProfileMatcherParams(const ProfileMatcherParams& p)
            : ignore_context(p.ignore_context),
              weight_center(p.weight_center),
              weight_decay(p.weight_decay)
    { }

    virtual ~ProfileMatcherParams()
    { }

    bool ignore_context;
    float weight_center;
    float weight_decay;
};



template< class Alphabet_T>
class ProfileMatcher
{
  public:
    // Constructs a profile matcher that scores only central column and ignores context.
    ProfileMatcher(int num_cols = 1);
    // Constructs a profile matcher with positional window weights.
    ProfileMatcher(int num_cols, const ProfileMatcherParams& params);
    // Don't forget to delete weights vector!
    ~ProfileMatcher();

    // Matches a context profile with a profile window centered at given index.
    double operator() (const ContextProfile<Alphabet_T>& profile,
                       const CountsProfile<Alphabet_T>& counts,
                       int index) const;
    // Matches a context profile with a sequence window centered at given index.
    double operator() (const ContextProfile<Alphabet_T>& profile,
                       const Sequence<Alphabet_T>& seq,
                       int index) const;
    // Calculates the effective number of columns for matching, that is the sum of weights.
    float num_eff_cols() const;
    // Initializes positional window weights.
    void init_weights(const ProfileMatcherParams& params);

    // TODO: Implement matching methods against counts profiles or sequences of fixed length W.

  private:
    // Disallow copy and assign
    ProfileMatcher(const ProfileMatcher&);
    void operator=(const ProfileMatcher&);

    // Number of columns in context profiles.
    int num_cols_;
    // Index of central column in context profiles.
    int center_;
    // Positional window weights
    float* weights_;
};



template< class Alphabet_T>
inline ProfileMatcher<Alphabet_T>::ProfileMatcher(int num_cols)
        : num_cols_(1),
          center_((num_cols - 1) / 2),
          weights_(NULL)
{ }

template< class Alphabet_T>
inline ProfileMatcher<Alphabet_T>::ProfileMatcher(int num_cols, const ProfileMatcherParams& params)
        : num_cols_(num_cols),
          center_((num_cols - 1) / 2),
          weights_(NULL)
{
    if (num_cols_ % 2 != 1)
        throw Exception("Number of columns for profile matching should be odd but is %i!", num_cols_);
    init_weights(params);
}

template< class Alphabet_T>
ProfileMatcher<Alphabet_T>::~ProfileMatcher()
{
    if (!weights_) delete weights_;
}

template< class Alphabet_T>
double ProfileMatcher<Alphabet_T>::operator() (const ContextProfile<Alphabet_T>& profile,
                                               const CountsProfile<Alphabet_T>& counts,
                                               int index) const
{
    double rv = 0.0f;
    if (weights_) {
        const int beg = std::max(0, index - center_);
        const int end = std::min(counts.num_cols() - 1, index + center_);
        for(int i = beg; i <= end; ++i) {
            int j = i - index + center_;
            double sum = 0.0f;
            for (int a = 0; a < profile.alphabet_size(); ++a)
                sum += counts[i][a] * profile[j][a];
            rv += weights_[j] * sum;
        }
    } else {
        for (int a = 0; a < profile.alphabet_size(); ++a)
            rv += counts[index][a] * profile[center_][a];
        }
    return pow(2.0, rv);
}

template< class Alphabet_T>
double ProfileMatcher<Alphabet_T>::operator() (const ContextProfile<Alphabet_T>& profile,
                                               const Sequence<Alphabet_T>& seq,
                                               int index) const
{
    double rv = 0.0f;
    if (weights_) {
        const int beg = std::max(0, index - center_);
        const int end = std::min(seq.length() - 1, index + center_);
        for(int i = beg; i <= end; ++i) {
            int j = i - index + center_;
            rv += weights_[j] * profile[j][seq[i]];
        }
    } else {
        rv = profile[center_][seq[index]];
    }
    return pow(2.0, rv);
}

template< class Alphabet_T>
void ProfileMatcher<Alphabet_T>::init_weights(const ProfileMatcherParams& params)
{
    if (params.ignore_context || num_cols_ == 1) return;

    if (!weights_) delete weights_;
    weights_ = new float[num_cols_];

    weights_[center_] = params.weight_center;
    for (int i = 1; i <= center_; ++i) {
        float weight = params.weight_center * pow(params.weight_decay, i);
        weights_[center_ - i] = weight;
        weights_[center_ + i] = weight;
    }
}

template< class Alphabet_T>
float ProfileMatcher<Alphabet_T>::num_eff_cols() const
{
    if (!weights_) return 1.0f;

    float sum = 0.0f;
    for (int i = 0; i < num_cols_; ++i) sum += weights_[i];
    return sum;
}

}  // cs

#endif
