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

template< class Alphabet_T>
class ProfileMatcher
{
  public:
    ProfileMatcher()
            : weights_(NULL),
              window_length_(1)
    { }

    ~ProfileMatcher()
    {
        if (!weights_) delete weights_;
    }

    void set_weights(int len, float weight_center, float weight_decay)
    {
        if (len % 2 != 1)
            throw Exception("Unable to set weights: window length should be odd but is %i!", len);

        if (!weights_) delete weights_;
        weights_ = new float[len];

        const int center = (len - 1) / 2;
        window_length_ = len;
        weights_[center] = weight_center;
        for (int i = 1; i <= center; ++i) {
            float weight = weight_center * pow(weight_decay, i);
            weights_[center - i] = weight;
            weights_[center + i] = weight;
        }
    }

    float num_eff_cols()
    {
        if (!weights_) return 1.0f;

        float sum = 0.0f;
        for (int i = 0; i < window_length_; ++i) sum += weights_[i];
        return sum;
    }

    double operator() (const ContextProfile<Alphabet_T>& profile,
                       const CountsProfile<Alphabet_T>& counts,
                       int index) const
    {
        double rv = 0.0f;
        const int center = profile.center();
        if (weights_) {
            const int beg = std::max(0, index - center);
            const int end = std::min(counts.num_cols() - 1, index + center);
            for(int i = beg; i <= end; ++i) {
                int j = i - index + center;
                double sum = 0.0f;
                for (int a = 0; a < profile.alphabet_size(); ++a)
                    sum += counts[i][a] * profile[j][a];
                rv += weights_[j] * sum;
            }
        } else {
            for (int a = 0; a < profile.alphabet_size(); ++a)
                rv += counts[index][a] * profile[center][a];
        }
        return pow(2.0, rv);
    }

    double operator() (const ContextProfile<Alphabet_T>& profile,
                       const Sequence<Alphabet_T>& seq,
                       int index) const
    {
        double rv = 0.0f;
        const int center = profile.center();
        if (weights_) {
            const int beg = std::max(0, index - center);
            const int end = std::min(seq.length() - 1, index + center);
            for(int i = beg; i <= end; ++i) {
                int j = i - index + center;
                rv += weights_[j] * profile[j][seq[i]];
            }
        } else {
            rv = profile[center][seq[index]];
        }
        return pow(2.0, rv);
    }

  private:
    // Disallow copy and assign
    ProfileMatcher(const ProfileMatcher&);
    void operator=(const ProfileMatcher&);

    // Positional window weights
    float* weights_;
    // Length of context window
    int window_length_;
};

}  // cs

#endif
