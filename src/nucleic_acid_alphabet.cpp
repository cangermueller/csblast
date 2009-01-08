/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleic_acid_alphabet.h"

namespace cs
{

NucleicAcidAlphabet::NucleicAcidAlphabet()
{ init(); }

NucleicAcidAlphabet::~NucleicAcidAlphabet()
{ }

NucleicAcidAlphabet* NucleicAcidAlphabet::instance()
{
    static NucleicAcidAlphabet inst;
    return &inst;
}

void NucleicAcidAlphabet::init_itoc()
{
    char itoc_arr[] = {'A','C','G','T','N'};
    itoc_.insert(itoc_.begin(), itoc_arr, itoc_arr + sizeof(itoc_arr)/sizeof(*itoc_arr));
}

}//cs
