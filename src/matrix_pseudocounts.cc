/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "matrix_pseudocounts.h"

#include "log.h"
#include "matrix.h"
#include "profile.h"
#include "counts_profile.h"
#include "shared_ptr.h"
#include "substitution_matrix.h"

namespace cs
{

MatrixPseudocounts::MatrixPseudocounts(const SubstitutionMatrix& m) : m_(m)
{}

void MatrixPseudocounts::add_to_sequence(const Sequence& seq, const AdmixtureCalculator& pca, Profile& p)
{
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to sequence ...";
    LOG(DEBUG1) << p;
    if (seq.alphabet() != p.alphabet())
        throw Exception("Cannot add substitution matrix pseudocounts: sequence and profile have different alphabets!");
    if (seq.length() != p.ncols())
        throw Exception("Cannot add substitution matrix pseudocounts: sequence and profile have different length!");

    const bool logspace = p.logspace();
    p.set_logspace(false);

    // add substitution matrix pseudocounts
    float tau = pca(1.0f);
    for(int i = 0; i < p.ncols(); ++i) {
        for(int a = 0; a < p.nalph(); ++a) {
            p[i][a] = (1.0f - tau) * (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * m_.r(a, seq[i]);
        }
    }

    normalize(p);
    if (logspace) p.transform_to_logspace();
    LOG(DEBUG1) << p;
}

void MatrixPseudocounts::add_to_profile(CountsProfile& p, const AdmixtureCalculator& pca)
{
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to profile ...";
    LOG(DEBUG1) << p;

    const bool logspace = p.logspace();
    if (logspace) p.transform_to_linspace();

    // copy original frequencies to matrix f
    matrix<float> f(p.ncols(), p.nalph(), 0.0f);
    for (int i = 0; i < p.ncols(); ++i)
        for(int a = 0; a < p.nalph(); ++a)
            f[i][a] = p[i][a];

    // add substitution matrix pseudocounts
    for(int i = 0; i < p.ncols(); ++i) {
        float tau = pca(p.neff(i));
        for(int a = 0; a < p.nalph(); ++a) {
            float sum = 0.0f;
            for(int b = 0; b < p.nalph(); ++b)
                sum += m_.r(a,b) * f[i][b];
            p[i][a] = (1.0f - tau) * f[i][a] + tau * sum;
        }
    }

    if (logspace) p.transform_to_logspace();
    normalize(p);
    LOG(DEBUG1) << p;
}


}  // cs
