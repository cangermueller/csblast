/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleic_acid_alphabet.h"

namespace cs
{

NucleicAcidAlphabet::NucleicAcidAlphabet()
{
    init();
}


NucleicAcidAlphabet::~NucleicAcidAlphabet()
{ }


NucleicAcidAlphabet* NucleicAcidAlphabet::instance()
{
    static NucleicAcidAlphabet inst;
    return &inst;
}


std::vector<char> NucleicAcidAlphabet::itoc() const
{
    char itoc_arr[] = {'A','C','G','T','N'};
    std::vector<char> itoc(itoc_arr, itoc_arr + sizeof(itoc_arr)/sizeof(*itoc_arr));
    return itoc;
}

}//cs
