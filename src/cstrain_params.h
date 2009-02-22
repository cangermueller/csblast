#ifndef CS_CSTRAIN_PARAMS_H
#define CS_CSTRAIN_PARAMS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Parameters for HMM training with Baum-Welch algorithm.

#include <iostream>
#include <string>

#include "exception.h"
#include "forward_backward_algorithm.h"
#include "getopt_pp.h"
#include "util.h"

using namespace GetOpt;

namespace cs
{

template<class Alphabet_T>
class CSTrainParams : public ForwardBackwardParams
{
  public:
    CSTrainParams();
    virtual ~CSTrainParams() { }

    int num_states() { return num_states_; }
    CSTrainParams& num_states(int n) { num_states_ = n; return *this; }

    void check();
    void parse(GetOpt_pp& options);
    std::ostream& usage(std::ostream& out);

  protected:
    int num_states_;
    int log_level_;
};



template<class Alphabet_T>
CSTrainParams<Alphabet_T>::CSTrainParams()
        : num_states_(0),
          log_level_(DEBUG)
{ }

template<class Alphabet_T>
std::ostream& CSTrainParams<Alphabet_T>::usage(std::ostream& out)
{
    out << "Train an HMM of context-states on a dataset of alignments.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain -i <alifile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << strprintf("  %-30s %s\n",          "-K, --num-states <int>", "Number of states in the HMM to be trained");
    out << strprintf("  %-30s %s (def=%i)\n", "-l, --log-max-level <int>", "Maximal reporting level for logging", log_level_);

    return out;
}

template<class Alphabet_T>
void CSTrainParams<Alphabet_T>::parse(GetOpt_pp& options)
{
    options >> Option('K', "num-states", num_states_, num_states_);
    options >> Option('l', "max-log-level", log_level_, log_level_);
    Log::reporting_level() = Log::from_integer(log_level_);
}

template<class Alphabet_T>
void CSTrainParams<Alphabet_T>::check()
{
    if (num_states() == 0) throw cs::Exception("No value for number of HMM states provided!");
}

}  // cs

#endif
