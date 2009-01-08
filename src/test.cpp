/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include <iostream>
#include <cstddef>
#include <vector>
#include <cmath>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"
#include "sequence_profile.h"

using std::cout;
using std::endl;

int main(int argc, const char **argv) {
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();

    for(cs::SequenceAlphabet::const_iterator iter=aa->begin(); iter != aa->end(); ++iter)
        cout << *iter << endl;
    cout << endl;
    for(cs::SequenceAlphabet::const_iterator iter=na->begin(); iter != na->end(); ++iter)
        cout << *iter << endl;

    cs::SequenceProfile profile(10, na);

    cout << endl << profile.ncols() << endl << profile.nalph() << endl;

    profile(0,4) = 0.5;
    cout << endl << profile(0,4) << endl << endl;

    for(cs::SequenceAlphabet::const_iterator iter=profile.alphabet().begin(); iter != profile.alphabet().end(); ++iter)
        cout << *iter << endl;

    return 0;
}
