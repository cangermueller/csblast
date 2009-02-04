/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "amino_acid.h"

#include "alphabet.h"

namespace cs
{

const char AminoAcid::amino_acids_[] = "ARNDCQEGHILKMFPSTWYV";

AminoAcid::AminoAcid() : Alphabet(20, 'X')
{
    init();
    // Conversion of non-standard amino acids to integers
    set_ctoi('O', any());
    set_ctoi('J', any());
    set_ctoi('U', ctoi('C'));
    set_ctoi('B', ctoi('D'));
    set_ctoi('Z', ctoi('E'));
}

AminoAcid& AminoAcid::instance()
{
    static AminoAcid inst;
    return inst;
}

};  // cs
