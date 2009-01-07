/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include <iostream>
#include <cstddef>
#include <vector>
#include <cmath>

#include "Amino_acid_alphabet.h"
#include "Nucleic_acid_alphabet.h"
#include "Sequence_profile.h"

using std::cout;
using std::endl;

int main(int argc, const char **argv) {
    Amino_acid_alphabet* aa = Amino_acid_alphabet::instance();
    Nucleic_acid_alphabet* na = Nucleic_acid_alphabet::instance();

    for(Sequence_alphabet::const_iterator iter=aa->begin(); iter != aa->end(); ++iter)
        cout << *iter << endl;
    cout << endl;
    for(Sequence_alphabet::const_iterator iter=na->begin(); iter != na->end(); ++iter)
        cout << *iter << endl;

    Sequence_profile profile(10, na);

    cout << endl << profile.ncols() << endl << profile.nalph() << endl;

    profile(0,4) = 0.5;
    cout << endl << profile(0,4) << endl << endl;

    for(Sequence_alphabet::const_iterator iter=profile.alphabet().begin(); iter != profile.alphabet().end(); ++iter)
        cout << *iter << endl;

    return 0;
}
