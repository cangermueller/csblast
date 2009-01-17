/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleotide_alphabet.h"

namespace cs
{

const char NucleotideAlphabet::nucleotides_[] = "ACGT";

NucleotideAlphabet::NucleotideAlphabet()
        : SequenceAlphabet(4, 'N', '-')
{ init(); }

NucleotideAlphabet* NucleotideAlphabet::instance()
{
    static NucleotideAlphabet inst;
    return &inst;
}

}//cs

