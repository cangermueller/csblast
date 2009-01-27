/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "profile.h"

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "exception.h"
#include "log.h"
#include "profile.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

Profile::Profile(const SequenceAlphabet* alphabet)
        : logspace_(false),
          alphabet_(alphabet)
{}

Profile::Profile(int ncols,
                 const SequenceAlphabet* alphabet)
        : data_(ncols, alphabet->size(), 0.0f),
          logspace_(false),
          alphabet_(alphabet)
{}

Profile::Profile(std::istream& in, const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ read(in); }

Profile::Profile(const Profile& other,
                 int index,
                 int length)
        : data_(length, other.nalph(), 0.0f),
          alphabet_(other.alphabet_)
{
    if (index + length > other.ncols())
        throw Exception("Arguments index=%i and length=%i of sub-profile are out of bounds!", index, length);
    for (int i = 0; i < ncols(); ++i)
        for (int a = 0; a < nalph(); ++a)
            data_[i][a] = other[i+index][a];
}

std::vector< shared_ptr<Profile> > Profile::readall(std::istream& in,
                                                 const SequenceAlphabet* alphabet)
{
    std::vector< shared_ptr<Profile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<Profile> p(new Profile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void Profile::transform_to_logspace()
{
    if (!logspace_) {
        for(int i = 0; i < ncols(); ++i)
            for(int a = 0; a < nalph(); ++a)
                data_[i][a] = data_[i][a] == 0.0f ? -std::numeric_limits<float>::infinity() : log2(data_[i][a]);
        logspace_ = true;
    }
}

void Profile::transform_to_linspace()
{
    if (logspace_) {
        for(int i = 0; i < ncols(); ++i)
            for(int a = 0; a < nalph(); ++a)
                data_[i][a] = pow(2.0, data_[i][a]);
        logspace_ = false;
    }
}

void Profile::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading profile from stream ...";

    // Check if stream actually contains a serialized profile
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find(class_identity()) == std::string::npos)
        throw Exception("Bad format: serialized profile does not start with '%s'!", class_identity().c_str());

    read_header(in);
    read_body(in);

    LOG(DEBUG1) << *this;
}

void Profile::read_header(std::istream& in)
{
    std::string tmp;

    // Read ncols
    int ncols = 0;
    if (getline(in, tmp) && tmp.find("ncols") != std::string::npos)
        ncols = atoi(tmp.c_str() + 5);
    else
        throw Exception("Bad format: serialized profile does not contain 'ncols' record!");

    // Read nalph
    int nalph = 0;
    if (getline(in, tmp) && tmp.find("nalph") != std::string::npos)
        nalph = atoi(tmp.c_str() + 5);
    else
        throw Exception("Bad format: serialized profile does not contain 'nalph' record!");
    if (nalph != alphabet_->size())
        throw Exception("Bad format: nalph=%i does not fit with provided alphabet!", nalph);

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    resize(ncols, nalph);
}

void Profile::read_body(std::istream& in)
{
    std::string tmp;
    std::vector<std::string> tokens;
    int i = 0;
    getline(in, tmp);  // skip alphabet description line
    while (getline(in, tmp)) {
        if (tmp.empty()) continue;
        if (tmp.length() > 1 && tmp[0] == '/' && tmp[1] == '/') break;

        split(tmp, '\t', tokens);
        i = atoi(tokens[0].c_str()) - 1;
        for (int a = 0; a < nalph(); ++a) {
            float log_p = tokens[a+1][0] == '*' ? std::numeric_limits<int>::max() : atoi(tokens[a+1].c_str());
            data_[i][a] = (logspace_ ? -log_p / kScaleFactor : pow(2.0, -log_p / kScaleFactor)) ;
        }
        tokens.clear();
    }
    if (i != ncols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, ncols());
}

void Profile::write(std::ostream& out) const
{
    out << class_identity() << std::endl;
    write_header(out);
    write_body(out);
}

void Profile::write_header(std::ostream& out) const
{
    // print dimensions
    out << "ncols\t\t" << ncols() << std::endl;
    out << "nalph\t\t" << nalph() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;
}

void Profile::write_body(std::ostream& out) const
{
    out << "\t";
    alphabet_->print(out, "\t");
    out << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < nalph(); ++a) {
            float logval = logspace_ ? data_[i][a] : log2(data_[i][a]);
            if (-logval == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * kScaleFactor);
        }
        out << std::endl;
    }
}

void Profile::print(std::ostream& out) const
{
    const int kWidth = 6;
    std::ios_base::fmtflags flags = out.flags(); // save old flags

    out << std::string(2*kWidth-2, ' ');
    alphabet_->print(out, std::string(kWidth-1, ' '));
    out << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << std::left << std::setw(kWidth-1) << i+1;
        for (int a = 0; a < nalph(); ++a)
            out << std::right << std::setw(kWidth) << std::fixed << std::setprecision(3)
                << (logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
        out << std::endl;
    }

    out.flags(flags);
}

void Profile::resize(int ncols, int nalph)
{
    if (ncols == 0 || nalph == 0)
        throw Exception("Bad dimensions for profile resizing: ncols=%i nalph=%i", ncols, nalph);
    data_.resize(ncols, nalph);
}

void reset(Profile& profile, float value)
{
    const int ncols = profile.ncols();
    const int nalph = profile.nalph();
    for(int i = 0; i < ncols; ++i)
        for(int a = 0; a < nalph; ++a)
            profile[i][a] = value;
}

void normalize(Profile& profile, float value)
{
    const bool logspace = profile.logspace();
    if (logspace) profile.transform_to_linspace();

    const int ncols = profile.ncols();
    const int nalph  = profile.nalph();
    for (int i = 0; i < ncols; ++i) {
        float sum = 0.0f;
        for (int a = 0; a < nalph; ++a) sum += profile[i][a];
        if (sum != 0.0f) {
            float fac = value / sum;
            for (int a = 0; a < nalph; ++a) profile[i][a] *= fac;
        } else {
            throw Exception("Unable to normalize profile to one. Sum of column %i is zero!", i);
        }
    }

    if (logspace) profile.transform_to_logspace();
}

}//cs


