// Copyright 2009, Andreas Biegert

#ifndef CS_LR_PSEUDOCOUNTS_INL_H_
#define CS_LR_PSEUDOCOUNTS_INL_H_

#include "lr_pseudocounts.h"


namespace cs {


template<class Abc>
void LrPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const {
	assert_eq(seq.length(), p.length());
  LOG(INFO) << "Adding lr pseudocounts to sequence ...";
	LrKmerBuilder<Abc> b(params_.cfg());
	size_t w = params_.cfg().wlen / 2;
	for (size_t i = 0; i < seq.length(); i++) {
		vector<long> kmers = b(seq, i - w);
		float e[Abc::kSize];
		float e_sum = 0;
		for (size_t a = 0; a < Abc::kSize; a++) {
			e[a] = exp(params_.sum(kmers, a));
			e_sum += e[a];
		}
		for (size_t a = 0; a < Abc::kSize; a++)
			p[i][a] = e[a] / e_sum;
	}
}

template<class Abc>
void LrPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const {

 	assert_eq(cp.counts.length(), p.length());
  LOG(INFO) << "Adding lr pseudocounts to profile ...";
  Sequence<Abc> seq = ConsensusSeq(cp);
	LrKmerBuilder<Abc> b(params_.cfg());
	size_t w = params_.cfg().wlen / 2;
	for (size_t i = 0; i < cp.counts.length(); i++) {
		vector<long> kmers = b(seq, i - w);
		float e[Abc::kSize];
		float e_sum = 0;
		for (size_t a = 0; a < Abc::kSize; a++) {
			e[a] = exp(params_.sum(kmers, a));
			e_sum += e[a];
		}
		for (size_t a = 0; a < Abc::kSize; a++)
			p[i][a] = e[a] / e_sum;
	}
}

template<>
Sequence<AA> LrPseudocounts<AA>::ConsensusSeq(const CountProfile<AA>& cp) const {
	return ConsensusSequence(cp, BlosumMatrix());
}

template<>
Sequence<Dna> LrPseudocounts<Dna>::ConsensusSeq(const CountProfile<Dna>& cp) const {
	return ConsensusSequence(cp, TamuraNeiMatrix());
}



}  // namespace cs

#endif  // CS_LR_PSEUDOCOUNTS_INL_H_
