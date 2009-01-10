/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence.h"

namespace cs
{

Sequence::Sequence(int length,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet),
          sequence_(length)
{}

Sequence::Sequence(const std::string& header,
                   const std::string& sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{
    //adjust capacity of header and sequence a minimum
    std::string().swap(header_);
    std::vector<char>().swap(sequence_);
    sequence_.reserve(sequence.size());

    header_.insert(header_.begin(), header.begin() + (header[0]=='>' ? 1 : 0), header.end());
    sequence_.insert(sequence_.begin(), sequence.begin(), sequence.end());

    check_and_convert();
}

Sequence::Sequence(std::istream& in,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ in >> *this; }

Sequence::~Sequence()
{}

std::vector<Sequence*> Sequence::read(std::istream& in,
                                      const SequenceAlphabet* alphabet)
{
    std::vector<Sequence*> sequences;
    while (in.good()) {
        Sequence* seq = new Sequence(in, alphabet);
        sequences.push_back(seq);
    }

    return sequences;
}

// Convert the sequence in character representation to integer representation. This also
// involves removing whitespace and checking for invalid characters. In the latter case
// an exception is thrown.
void Sequence::check_and_convert()
{
    //strip whitespace and newlines from sequence.
    sequence_.erase(remove_if(sequence_.begin(), sequence_.end(), isspace),
                    sequence_.end());
    //validate each character and convert to integer representation
    int len = sequence_.size();
    for (int i = 0; i < len; ++i)
        if (alphabet_->valid(sequence_[i])) {
            sequence_[i] = alphabet_->ctoi(sequence_[i]);
        } else {
            char c = sequence_[i];
            sequence_.clear();
            throw MyException("Invalid character %c at position %i of sequence '%s'", c, i, header_.c_str());
        }
}

// Initializes the sequence object from given stream with a sequence in FASTA format.
// Note that the sequence undergoes a validation check: White space is silently
// removed whereas invalid characters (according to the alphabet) trigger an
// execption. After the initialization the stream points to begin of next FASTA
// sequence or EOF.
void Sequence::init(std::istream& in)
{
    //clear and minimize capacity of header and sequence
    std::string().swap(header_);
    std::vector<char>().swap(sequence_);

    const int kBufferSize = 2097152; //2MB
    std::string buffer;
    buffer.reserve(kBufferSize);
    //read header
    if (getline(in, buffer)) {
        if (buffer.empty() ||  buffer[0] != '>')
            throw MyException("Bad format: first line of FASTA sequence does not start with '>' character!");
        header_.insert(header_.begin(), buffer.begin() + 1, buffer.end());
    } else {
        throw MyException("Failed to read from FASTA formatted input stream!");
    }
    //read sequence
    while(in.peek() != '>' && getline(in, buffer))
        sequence_.insert(sequence_.end(), buffer.begin(), buffer.end());

    check_and_convert();
}

std::istream& operator>> (std::istream& in, Sequence& sequence)
{
    sequence.init(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const Sequence& sequence)
{
    const int kLineLength = 80;

    out << '>' << sequence.header() << std::endl;
    int len = sequence.length();
    for (int i = 0; i < len; ++i) {
        out << sequence.chr(i);
        if ((i+1) % kLineLength == 0) out << std::endl;
    }

    return out;
}

}//cs
