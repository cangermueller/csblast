// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-24		

#ifndef CS_LR_PARAMS_INL_H_
#define CS_LR_PARAMS_INL_H_

#include "lr_params.h"


namespace cs {


// LrParams


template<class Abc>
LrCfg<Abc> LrParams<Abc>::ReadCfg(FILE *fin) {
    if (!StreamStartsWith(fin, "LR"))
        throw Exception("Stream does not start with class id 'LR'!");

	size_t wlen = 0;
	size_t klen = 0;
    char buffer[KB];
    if (fgetline(buffer, KB, fin)) 
        wlen = static_cast<size_t>(ReadInt(buffer, "WLEN", "Unable to parse LR 'WLEN'!"));
    if (fgetline(buffer, KB, fin))
        klen = static_cast<size_t>(ReadInt(buffer, "KLEN", "Unable to parse LR 'KLEN'!"));
	if (wlen == 0 || klen == 0) 
		throw Exception("Stream has no valid window configuration!");
	return LrCfg<Abc>(wlen, klen);
}

template<class Abc>
void LrParams<Abc>::Read(FILE* fin) {
    // Parse and check header information
    if (!StreamStartsWith(fin, "LR"))
        throw Exception("Stream does not start with class id 'LR'!");

    char buffer[KB];
    if (fgetline(buffer, KB, fin)) {
        if (static_cast<unsigned int>(ReadInt(buffer, "WLEN", "Unable to parse LR 'WLEN'!")) != cfg_.wlen) 
			throw Exception("WLEN must be %i!", cfg_.wlen);
	}
    if (fgetline(buffer, KB, fin)) {
        if (static_cast<unsigned int>(ReadInt(buffer, "KLEN", "Unable to parse LR 'LEN'!")) != cfg_.klen)
			throw Exception("KLEN must be %i!", cfg_.klen);
	}
    if (fgetline(buffer, KB, fin)) {
        if (static_cast<unsigned int>(ReadInt(buffer, "ALPH", "Unable to parse LR 'ALPH'!")) != Abc::kSize) 
			throw Exception("ALPH must have %i letters!", Abc::kSize);
    }

    // Read params
    fgetline(buffer, KB, fin);
    for (size_t p = 0; p < cfg_.nparams_let; p++) {
    	fgetline(buffer, KB, fin);
		const char* ptr = buffer;
		for (size_t a = 0; a < Abc::kSize; a++) {
			int i = strastoi(ptr);
			if (i == INT_MIN) 
				throw Exception("Missing parameters in LR file!");
			(*params_)[a][p] = static_cast<float>(i) / kPrecision;
		}
	}
}

template<class Abc>
void LrParams<Abc>::Write(FILE* fout) const {
    // Write header
    fputs("LR\n", fout);
    fprintf(fout, "WLEN\t%3d\n", static_cast<int>(cfg_.wlen));
    fprintf(fout, "KLEN\t%3d\n", static_cast<int>(cfg_.klen));
    fprintf(fout, "ALPH\t%3d\n",  static_cast<int>(Abc::kSize));

    // Print letters
    for (size_t a = 0; a < Abc::kSize; a++) 
		fprintf(fout, "%10c", Abc::kIntToChar[a]);
    fprintf(fout, "\n");
    // Serialize states
	for (size_t p = 0; p < cfg_.nparams_let; p++) {
		for (size_t a = 0; a < Abc::kSize; a++) 
			fprintf(fout, "%10i", static_cast<int>((*params_)[a][p] * kPrecision));
    	fprintf(fout, "\n");
    }

}

template<class Abc>
void LrParams<Abc>::operator= (const LrParams<Abc>& p) {
	if (cfg_ != p.cfg())
		throw Exception("wlen or klen differ!");
	*params_ = p.params(); 
}

template<class Abc>
inline float LrParams<Abc>::sum(std::vector<long> kmers, size_t a) const {
	float sum = 0;
	std::vector<long>::iterator it = kmers.begin();
	while (it != kmers.end()) {
		sum += (*params_)[a][*it];
		it++;
	}
	return sum;
}


// GaussianLrParamsInit


template<class Abc>
void GaussianLrParamsInit<Abc>::operator() (LrParams<Abc>& params) const {
    Gaussian gauss(0, sigma_, seed_);
	float* p = params.params().begin();
	for (size_t i = 0; i < params.cfg().nparams; i++, p++) 
		*p = gauss();
}


}  // namespace cs

#endif  // CS_LR_PARAMS_INL_H_
