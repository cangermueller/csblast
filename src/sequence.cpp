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
                   const std::vector<char>& sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ init(header, sequence); }

Sequence::Sequence(std::istream& in,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ in >> *this; }

Sequence::~Sequence()
{}

// Initializes sequence object with given header string and sequence vector.
// Note that the sequence undergoes a validation check: White space is silently
// removed whereas invalid characters (according to the alphabet) trigger an
// execption.
void Sequence::init(const std::string& header, const std::vector<char>& sequence)
{
    std::string().swap(header_); //clear and minimize capacity
    header_.insert(header_.begin(), header.begin() + (header[0]=='>' ? 1 : 0), header.end());

    std::vector<char>().swap(sequence_); //clear and minimize capacity
    sequence_.reserve(sequence.size());
    sequence_.insert(sequence_.begin(), sequence.begin(), sequence.end());
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

std::istream& operator>> (std::istream& in, Sequence& sequence)
{
    std::string data((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());

    size_t i = data.find_first_of('>');
    if (i == std::string::npos) throw MyException("Bad format: FASTA input does not contain '>'!");
    size_t j = data.find_first_of('\n', i);
    if (j == std::string::npos) throw MyException("Bad format: FASTA header does not terminate with newline!");
    size_t k = data.find_first_of('>', j); //check if there is another FASTA sequence after the first

    std::string header(data.begin()+i, data.begin()+j);
    std::vector<char> raw_sequence(data.begin()+j, k == std::string::npos ? data.end() : data.begin()+k);
    sequence.init(header, raw_sequence);

    return in;
}

std::ostream& operator<< (std::ostream& out, const Sequence& sequence)
{
    const int kLineLength = 80;

    out << '>' << sequence.header() << std::endl;
    int len = sequence.length();
    for (int i = 0; i < len; ++i) {
        out << sequence.chr(i);
        if (i+1 % kLineLength == 0) out << std::endl;
    }

    return out;
}

}//cs
