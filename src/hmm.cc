/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "hmm.h"

#include <algorithm>
#include <functional>
#include <list>
#include <vector>
#include <utility>

#include "exception.h"
//#include "initializer.h"
#include "log.h"
#include "profile.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"

namespace cs
{

const char HMM::kClassIdentity[] = "HMM";

void HMM::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading HMM from stream ...";

    // Check if stream actually contains a serialized HMM
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find(kClassIdentity) == std::string::npos)
        throw Exception("Bad format: serialized HMM does not start with '%s' keyword!", kClassIdentity);

    // Read number of states
    int size = 0;
    if (getline(in, tmp) && tmp.find("size") != std::string::npos)
        size = atoi(tmp.c_str() + 4);
    else
        throw Exception("Bad format: serialized profile does not contain 'size' record!");

    states_.clear();
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<State> state(new State(in, alphabet_));
        states_.push_back(state);
    }

    LOG(DEBUG1) << *this;
}

void HMM::write(std::ostream& out) const
{
    out << kClassIdentity << std::endl;
    out << "size\t" << size() << std::endl;
    for (std::pair<const_state_iterator, const_state_iterator> sp = states(); sp.first != sp.second; ++sp.first)
        (**sp.first).write(out);
}

void HMM::sparsify(float threshold)
{
    std::pair<state_iterator, state_iterator> sp = states();
    for (std::pair<state_iterator, state_iterator> sp = states(); sp.first != sp.second; ++sp.first) {
        (*sp.first)->in_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
        (*sp.first)->out_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
    }
}



void HMM::State::read_header(std::istream& in)
{
    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("index") != std::string::npos)
        index_ = atoi(tmp.c_str() + 5);

    Profile::read_header(in);
}

void HMM::State::read_body(std::istream& in)
{
    // TODO
}

void HMM::State::write_header(std::ostream& out) const
{
    out << "index\t" << index_ << std::endl;
    Profile::write_header(out);
}

void HMM::State::write_body(std::ostream& out) const
{
    // TODO
}

void HMM::State::print(std::ostream& out) const
{
    out << "State " << index_ << ":" << std::endl;
    Profile::print(out);
}

}  // cs
