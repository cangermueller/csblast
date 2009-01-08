/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence.h"

namespace cs
{

Sequence::Sequence(int length,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet), sequence_(length)
{}

Sequence::Sequence(const std::string header,
                   const std::string sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ init(header.c_str(), sequence.c_str()); }

Sequence::Sequence(const char* header,
                   const char* sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ init(header, sequence); }

Sequence::~Sequence()
{}

// TODO: implement Sequence::init
void Sequence::init(const char* header, const char* sequence)
{ return; }

// TODO: implement input operator >> for class Sequence
std::istream& operator>> (std::istream& i, Sequence& sequence)
{ return i; }

// TODO: implement output operator << for class SequenceProfile
std::ostream& operator<< (std::ostream& o, const Sequence& sequence)
{ return o; }

}//cs
