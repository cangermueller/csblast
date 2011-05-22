// Copyright 2009, Andreas Biegert

#ifndef CS_CRF_INL_H_
#define CS_CRF_INL_H_

#include "crf.h"

#include "crf_state-inl.h"
#include "pseudocounts-inl.h"

namespace cs {

template<class Abc>
Crf<Abc>::Crf(size_t size, size_t wlen)
        : wlen_(wlen), states_(size, CrfState<Abc>(wlen)) {}

template<class Abc>
Crf<Abc>::Crf(FILE* fin)
        : wlen_(0), states_() {
    Read(fin);
}

template<class Abc>
Crf<Abc>::Crf(size_t size, size_t wlen, const CrfInit<Abc>& init)
        : wlen_(wlen), states_(size, CrfState<Abc>(wlen)) {
    init(*this);
}

template<class Abc>
inline void Crf<Abc>::SetState(size_t k, const CrfState<Abc>& p) {
    assert_eq(wlen(), p.context_weights.length());
    assert(k < size());
    states_[k] = p;
}

template<class Abc>
void Crf<Abc>::Read(FILE* fin) {
    // Parse and check header information
    if (!StreamStartsWith(fin, "CRF"))
        throw Exception("Stream does not start with class id 'CRF'!");

    char buffer[KB];
    size_t size = 0;
    if (fgetline(buffer, KB, fin))
        size = ReadInt(buffer, "SIZE", "Unable to parse CRF 'SIZE'!");
    if (fgetline(buffer, KB, fin))
        wlen_ = ReadInt(buffer, "LENG", "Unable to parse CRF 'LENG'!");

    // Read context states
    states_.Resize(size);
    size_t k = 0;
    for (; k < size && !feof(fin); ++k) {
        states_[k] = CrfState<Abc>(fin);
    }
    LOG(DEBUG1) << *this;
    if (k != size)
        throw Exception("Serialized CRF should have %i states but actually has %i!",
                        size, k);
}

template<class Abc>
void Crf<Abc>::Write(FILE* fout) const {
    // Write header
    fputs("CRF\n", fout);
    fprintf(fout, "SIZE\t%d\n", static_cast<int>(size()));
    fprintf(fout, "LENG\t%d\n", static_cast<int>(wlen()));
    // Serialize states
    for (size_t k = 0; k < states_.size(); ++k) states_[k].Write(fout);
}


template<class Abc, class TrainingPair>
void SamplingCrfInit<Abc, TrainingPair>::operator() (Crf<Abc>& crf) const {
    LOG(DEBUG) << "Initializing CRF with by sampling " << crf.size()
               << " profile windows from training set ...";

    if (trainset_.size() < crf.size())
      throw Exception("CRF initialization requires %zu training pairs but the"
        "training set only contains %zu training pairs!", crf.size(), trainset_.size());
          
    Ran ran(seed_);
    Vector<bool> used(trainset_.size(), false);

    size_t k = 0;
    while (k < crf.size()) {
        size_t r = ran(crf.size());
        if (!used[r]) {
          if (crf.wlen() != trainset_[r].x.length())
            throw Exception("CRF window length does not match the training set window length!");
          Profile<Abc> prof = pc_.AddTo(trainset_[r].x, admix_);
          // Avoid null columns due to ANY char. TODO: Delete this block if ANY is serialized
          for (size_t i = 0; i < prof.length(); ++i) {
            bool null = true;
            for (size_t a = 0; a < Abc::kSize; ++a) {
              if (prof[i][a] > 0) {
                null = false;
                break;
              }
            }
            if (null) {
              double tau = admix_(1.0);
              for (size_t a = 0; a < Abc::kSize; ++a)
                prof[i][a] = tau * sm_.r(a, Abc::kAny);
            }
          }
          // Initialize CrfSate so that pc_weights are equal to context_weights[center]
          CrfState<Abc> s(1.0 / crf.size(), prof, trainset_[r].y, 
              weight_center_, weight_decay_); 
          const size_t c = crf.center();
          for (size_t a = 0; a < Abc::kSize; ++a)
            s.pc_weights[a] = s.context_weights[c][a];            
          UpdatePseudocounts(s);
              
          crf.SetState(k++, s);
          used[r] = true;
        }
    }

    LOG(DEBUG) << crf;
}


template<class Abc>
void GaussianCrfInit<Abc>::operator() (Crf<Abc>& crf) const {
    Gaussian gauss(0, sigma_, seed_);

    for (size_t k = 0; k < crf.size(); ++k) {
        CrfState<Abc> s(crf.wlen());
        s.bias_weight = gauss();
        for (size_t j = 0; j < crf.wlen(); ++j) {
            for (size_t a = 0; a < Abc::kSize; ++a)
                s.context_weights[j][a] = gauss();
            s.context_weights[j][Abc::kAny] = 0.0;
        }
        for (size_t a = 0; a < Abc::kSize; ++a)
            s.pc_weights[a] = log(sm_.p(a)) + gauss();
        s.pc_weights[Abc::kAny] = 0.0;
        UpdatePseudocounts(s);
        crf.SetState(k, s);
    }
}

}  // namespace cs

#endif  // CS_CRF_INL_H_
