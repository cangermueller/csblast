/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence.h"

#include <cctype>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "exception.h"
#include "log.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"

namespace cs
{

Sequence::Sequence(int length, const SequenceAlphabet* alphabet)
        : alphabet_(alphabet),
          seq_(length)
{}

Sequence::Sequence(std::istream& in, const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ read(in); }

Sequence::Sequence(const std::string& header,
                   const std::string& sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ init(header, sequence); }

std::vector< shared_ptr<Sequence> > Sequence::readall(std::istream& in,
                                                      const SequenceAlphabet* alphabet)
{
    std::vector< shared_ptr<Sequence> > sequences;
    while (in.good()) {
        shared_ptr<Sequence> p(new Sequence(in, alphabet));
        sequences.push_back(p);
    }

    return sequences;
}

void Sequence::init(std::string header, std::string sequence)
{
    // init header
    std::string().swap(header_);
    header_.append(header.begin() + (header[0]=='>' ? 1 : 0), header.end());

    //strip whitespace and newlines from sequence.
    sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace), sequence.end());
    //validate each character and convert to integer representation
    const int seqlen = sequence.length();
    seq_.resize(seqlen);
    for (int i = 0; i < seqlen; ++i) {
        char c = sequence[i];
        if (alphabet_->valid(c)) {
            seq_[i] = alphabet_->ctoi(c);
        } else {
            throw Exception("Invalid character %c at position %i of sequence '%s'", c, i, header_.c_str());
        }
    }
}

void Sequence::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading sequence from stream ...";
    std::string tmp;
    std::string header;
    std::string sequence;

    //read header
    if (getline(in, tmp)) {
        if (tmp.empty() ||  tmp[0] != '>')
            throw Exception("Bad format: first line of FASTA sequence does not start with '>' character!");
        header.append(tmp.begin() + 1, tmp.end());
    } else {
        throw Exception("Failed to read from FASTA formatted input stream!");
    }
    //read sequence
    while (in.peek() != '>' && getline(in, tmp))
        sequence.append(tmp.begin(), tmp.end());

    init(header, sequence);
    LOG(DEBUG1) << *this;
}

void Sequence::write(std::ostream& out, int width) const
{
    out << '>' << header_ << std::endl;
    for (int i = 0; i < length(); ++i) {
        out << chr(i);
        if ((i+1) % width == 0) out << std::endl;
    }
    if (length() % width != 0) out << std::endl;
}

}//cs
