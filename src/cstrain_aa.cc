/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include <iostream>

#include "amino_acid.h"
#include "cstrain_params.h"
#include "log.h"

int main(int argc, char* argv[])
{
    cs::CSTrainParams<cs::AminoAcid> params;
    GetOpt_pp options(argc, argv);

    if(argc < 2 || options >> OptionPresent('h', "help")) {
        params.usage(std::cout);
        return 2;
    }

    try {
        params.parse(options);
        params.check();

        std::cout << "num_states=" << params.num_states() << std::endl;
        std::cout << "log_level=" << Log::reporting_level() << std::endl;

    } catch(const std::exception& e) {
        LOG(ERROR) << e.what();
        std::cout << e.what() << std::endl;
        return 3;
    }
    return 0;
}
