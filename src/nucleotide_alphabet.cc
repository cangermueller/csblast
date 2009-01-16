/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleotide_alphabet.h"

namespace cs
{

const char NucleotideAlphabet::nucleotides_[] = "ACGTN";

NucleotideAlphabet::NucleotideAlphabet()
        : SequenceAlphabet(5)
{ init(); }

NucleotideAlphabet* NucleotideAlphabet::instance()
{
    static NucleotideAlphabet inst;
    return &inst;
}

}//cs

