// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-24		
//test 
#ifndef CS_LR_PARAMS_H_
#define CS_LR_PARAMS_H_


namespace cs {


// LrCfg


// Stores window configuration wlen and klen and the 
// number of parameters determined by wlen and klen.
template<class Abc>
struct LrCfg {
	LrCfg(size_t wlen_, size_t klen_) : 
		wlen(wlen_), 
		klen(klen_),
		nkmers(wlen - klen + 1),
		nparams_let(nkmers * pow(Abc::kSize, klen)),
		nparams(nparams_let * Abc::kSize) {}

	bool operator== (const LrCfg<Abc>& cfg) const { return wlen == cfg.wlen && klen == cfg.klen; }
	bool operator!= (const LrCfg<Abc>& cfg) const { return wlen != cfg.wlen || klen != cfg.klen; }

	const size_t wlen;			// window length
	const size_t klen;			// kmer length
	const size_t nkmers;		// number of kmers in a window: wlen - klen + 1
	const size_t nparams_let;	// number of kmer weights for a letter: |wa|
	const size_t nparams;		// total number of kmer weigths: sum |wa|
};


// LrParams


template<class Abc>
class LrParams;

// Strategy class for initializing lr params.
template<class Abc>
class LrParamsInit {
  public:
    LrParamsInit() {}
    virtual ~LrParamsInit() {}
    virtual void operator() (LrParams<Abc>& lrParams) const = 0;
};

// Stores kmer weigths w.
template<class Abc>
class LrParams {
  public:
    LrParams(const LrCfg<Abc>& cfg) : 
		cfg_(cfg),
		params_(new Matrix<float>(Abc::kSize, cfg_.nparams_let)) {}

    LrParams(const LrCfg<Abc>& cfg, const LrParamsInit<Abc>& init) :
		cfg_(cfg),
		params_(new Matrix<float>(Abc::kSize, cfg_.nparams_let)) {
	    init(*this); 
    }

    LrParams(FILE* fin) :
		cfg_(ReadCfg(fin)),
		params_(new Matrix<float>(Abc::kSize, cfg_.nparams_let)) {
		rewind(fin);
		Read(fin);
	}

	LrParams(const LrParams<Abc>& p) :
		cfg_(p.cfg()),
		params_(new Matrix<float>(p.params())) {}
 
    virtual ~LrParams() { delete params_; }

	// Access operator for wa.
    float* operator[] (int r) const { return (*params_)[r]; }

	// Assignment operator.
	void operator= (const LrParams<Abc>& p);

	// Reads parameters from a file with the same window config.
    void Read(FILE* fin);

	// Writes parameters to a file.
    void Write(FILE* fout) const;

	// Calculates sum of parameters for the given kmers and letter.
	inline float sum(std::vector<long> kmers, size_t a) const;

	// Returns window config.
	LrCfg<Abc> cfg() const { return cfg_; }

	// Returns parameter matrix.
	Matrix<float>& params() const { return *params_; }


  private:
	// Reads window config from file
	static LrCfg<Abc> ReadCfg(FILE *fin);

	static const size_t kPrecision = 10000000;	// Precision for serialization
	const LrCfg<Abc> cfg_;						// window config
  	Matrix<float>* const params_;				// parameters

};	// class LrParams


// GaussianLrParamsInit


// Strategy that initializes lr params by sammpling from gaussian distribution
template<class Abc>
class GaussianLrParamsInit : public LrParamsInit<Abc> {
  public:
    GaussianLrParamsInit(double sigma, unsigned int seed = 0) : 
		sigma_(sigma), seed_(seed) {}

    virtual ~GaussianLrParamsInit() {}

    virtual void operator() (LrParams<Abc>& LrParams) const;

  protected:
    double sigma_;
    unsigned int seed_;
};


}  // namespace cs

#endif	// CS_LR_PARAMS_H_
