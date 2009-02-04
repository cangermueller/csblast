/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleotide.h"

#include "alphabet.h"

namespace cs
{

const char Nucleotide::nucleotides_[] = "ACGT";

Nucleotide::Nucleotide() : Alphabet(4, 'N')
{
    init();
}

const Nucleotide& Nucleotide::instance()
{
    static Nucleotide inst;
    return inst;
}

};  // cs
