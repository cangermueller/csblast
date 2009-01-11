/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_profile.h"

namespace cs
{

SequenceProfile::SequenceProfile(int ncols, int ndim)
        : Profile(ncols, ndim),
          alphabet_(0)
{}

SequenceProfile::SequenceProfile(int ncols,
                                 const SequenceAlphabet* alphabet)
        : Profile(ncols, alphabet->size()-1),
          alphabet_(alphabet)
{}

SequenceProfile::SequenceProfile(const Sequence& sequence)
        : Profile(sequence.length(), sequence.alphabet().size()-1),
          alphabet_(&sequence.alphabet())
{
    const int cols = ncols();
    const int dim  = ndim();
    for(int i=0; i<cols; ++i) {
        for(int j=0; j<dim; ++j) (*this)(i,j) = 0.0f;
        (*this)(i, sequence(i));
    }
}

SequenceProfile::SequenceProfile(std::istream& in,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ in >> *this; }

SequenceProfile::~SequenceProfile()
{}

std::vector< SmartPtr<SequenceProfile> > SequenceProfile::read(std::istream& in,
                                                               const SequenceAlphabet* alphabet)
{
    std::vector< SmartPtr<SequenceProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        SmartPtr<SequenceProfile> p(new SequenceProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

SequenceProfile::SequenceProfile(const SequenceProfile& other,
                                 int index,
                                 int length)
        : Profile(length, other.ndim()),
          alphabet_(other.alphabet_)
{
    if (index + length > other.ncols())
        throw MyException("Arguments index=%i and length=%i for construction of sub-profile are out of bounds!", index, length);
    const int cols = ncols();
    const int dim  = ndim();
    for(int i=0; i<cols; ++i)
        for(int j=0; j<dim; ++j)
            (*this)(i,j) = other(i+index,j);
}

void SequenceProfile::init(std::istream& in)
{
    const int kBufferSize = 1048576; //1MB
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    //read column data records line by line
    std::vector< std::vector<int> > data;
    int i_prev = 0; //index of previous column
    while (in.getline(buffer, kBufferSize)) {
        int len = strlen(buffer);
        if (!strscn(buffer) || len > 0 && buffer[0] == '#') continue;
        if (len > 1 && buffer[0] == '/' && buffer[1] == '/') break;

        const char* ptr = buffer;
        std::vector<int> column;
        column.reserve(alphabet_->size()-1);
        int i = strtoi(ptr);
        if(i != i_prev+1)
            throw MyException("Bad format: column record %i is followed by column record %i!", i_prev, i);
        while ((ptr = strscn(ptr))) //read column values one by one
            column.push_back(strtoi_asterix(ptr));
        if(static_cast<int>(column.size()) != alphabet_->size()-1)
            throw MyException("Bad format: column record contains %i values but should have %i!",
                              static_cast<int>(column.size()), alphabet_->size()-1);
        data.push_back(column);
        i_prev = i;
    }
    const int nrows = data.size();
    const int ncols = alphabet_->size()-1;
    resize(nrows, ncols);

    //convert data to lin space anf fill matrix
    for (int i = 0; i < nrows; ++i)
        for (int j = 0; j < ncols; ++j)
            (*this)(i,j) = pow(2.0, static_cast<float>(-data[i][j]) / kScaleFactor);
}

std::istream& operator>> (std::istream& in, SequenceProfile& profile)
{
    profile.init(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile)
{
    const int ncols = profile.ncols();
    const int ndim  = profile.ndim();

    // print header with character alphabet
    out << "#";
    for (int j = 0; j < ndim; ++j)
        out << "\t" << profile.alphabet().itoc(j);
    out << std::endl;
    // print profile values in log representation
    for (int i = 0; i < ncols; ++i) {
        out << i+1;
        for (int j = 0; j < ndim; ++j) {
            double logval = log2(profile(i,j));
            if (-logval == std::numeric_limits<double>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * profile.kScaleFactor);
        }
        out << std::endl;
    }
    out << "//\n";

    return out;
}

}//cs


