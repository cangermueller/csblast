/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "blosum_matrix.h"

#include "amino_acid_alphabet.h"
#include "my_exception.h"
#include "substitution_matrix.h"

namespace
{

const float g_blosum45[] = {
    //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
    0.0181,
    0.0029,0.0130,
    0.0026,0.0020,0.0079,
    0.0027,0.0020,0.0031,0.0132,
    0.0015,0.0006,0.0007,0.0006,0.0094,
    0.0024,0.0025,0.0016,0.0017,0.0004,0.0057,
    0.0038,0.0031,0.0022,0.0047,0.0007,0.0032,0.0131,
    0.0063,0.0022,0.0033,0.0028,0.0010,0.0019,0.0025,0.0285,
    0.0013,0.0013,0.0013,0.0012,0.0003,0.0010,0.0014,0.0012,0.0058,
    0.0036,0.0016,0.0015,0.0014,0.0008,0.0013,0.0017,0.0019,0.0007,0.0124,
    0.0052,0.0029,0.0019,0.0022,0.0015,0.0022,0.0031,0.0032,0.0016,0.0093,0.0263,
    0.0037,0.0061,0.0028,0.0028,0.0008,0.0029,0.0045,0.0030,0.0013,0.0020,0.0031,0.0120,
    0.0016,0.0010,0.0007,0.0006,0.0004,0.0008,0.0009,0.0011,0.0006,0.0023,0.0041,0.0011,0.0026,
    0.0021,0.0014,0.0011,0.0010,0.0007,0.0007,0.0013,0.0017,0.0008,0.0031,0.0057,0.0015,0.0012,0.0124,
    0.0024,0.0013,0.0011,0.0015,0.0004,0.0011,0.0023,0.0022,0.0007,0.0016,0.0019,0.0020,0.0007,0.0009,0.0160,
    0.0062,0.0026,0.0031,0.0028,0.0011,0.0024,0.0032,0.0049,0.0013,0.0023,0.0032,0.0033,0.0010,0.0017,0.0020,0.0104,
    0.0041,0.0020,0.0025,0.0023,0.0010,0.0015,0.0025,0.0027,0.0009,0.0028,0.0038,0.0028,0.0011,0.0017,0.0019,0.0047,0.0086,
    0.0006,0.0004,0.0002,0.0002,0.0001,0.0003,0.0004,0.0006,0.0001,0.0005,0.0008,0.0005,0.0002,0.0008,0.0003,0.0004,0.0003,0.0053,
    0.0017,0.0014,0.0009,0.0011,0.0004,0.0010,0.0012,0.0014,0.0012,0.0019,0.0031,0.0015,0.0009,0.0034,0.0007,0.0015,0.0013,0.0008,0.0066,
    0.0055,0.0021,0.0016,0.0017,0.0012,0.0014,0.0023,0.0025,0.0008,0.0094,0.0088,0.0025,0.0022,0.0031,0.0016,0.0031,0.0038,0.0004,0.0019,0.0141 };

const float g_blosum62[] = {
    //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
    0.0215,
    0.0023,0.0178,
    0.0019,0.0020,0.0141,
    0.0022,0.0016,0.0037,0.0213,
    0.0016,0.0004,0.0004,0.0004,0.0119,
    0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,
    0.0030,0.0027,0.0022,0.0049,0.0004,0.0035,0.0161,
    0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,
    0.0011,0.0012,0.0014,0.0010,0.0002,0.0010,0.0014,0.0010,0.0093,
    0.0032,0.0012,0.0010,0.0012,0.0011,0.0009,0.0012,0.0014,0.0006,0.0184,
    0.0044,0.0024,0.0014,0.0015,0.0016,0.0016,0.0020,0.0021,0.0010,0.0114,0.0371,
    0.0033,0.0062,0.0024,0.0024,0.0005,0.0031,0.0041,0.0025,0.0012,0.0016,0.0025,0.0161,
    0.0013,0.0008,0.0005,0.0005,0.0004,0.0007,0.0007,0.0007,0.0004,0.0025,0.0049,0.0009,0.0040,
    0.0016,0.0009,0.0008,0.0008,0.0005,0.0005,0.0009,0.0012,0.0008,0.0030,0.0054,0.0009,0.0012,0.0183,
    0.0022,0.0010,0.0009,0.0012,0.0004,0.0008,0.0014,0.0014,0.0005,0.0010,0.0014,0.0016,0.0004,0.0005,0.0191,
    0.0063,0.0023,0.0031,0.0028,0.0010,0.0019,0.0030,0.0038,0.0011,0.0017,0.0024,0.0031,0.0009,0.0012,0.0017,0.0126,
    0.0037,0.0018,0.0022,0.0019,0.0009,0.0014,0.0020,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0014,0.0047,0.0125,
    0.0004,0.0003,0.0002,0.0002,0.0001,0.0002,0.0003,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0008,0.0001,0.0003,0.0003,0.0065,
    0.0013,0.0009,0.0007,0.0006,0.0003,0.0007,0.0009,0.0008,0.0015,0.0014,0.0022,0.0010,0.0006,0.0042,0.0005,0.0010,0.0009,0.0009,0.0102,
    0.0051,0.0016,0.0012,0.0013,0.0014,0.0012,0.0017,0.0018,0.0006,0.0120,0.0095,0.0019,0.0023,0.0026,0.0012,0.0024,0.0036,0.0004,0.0015,0.0196 };

const float g_blosum80[] = {
    //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
    0.0252,
    0.0020,0.0210,
    0.0016,0.0017,0.0166,
    0.0018,0.0013,0.0037,0.0262,
    0.0015,0.0003,0.0004,0.0003,0.0172,
    0.0017,0.0024,0.0014,0.0014,0.0003,0.0094,
    0.0028,0.0023,0.0019,0.0048,0.0003,0.0035,0.0208,
    0.0053,0.0015,0.0025,0.0022,0.0006,0.0011,0.0017,0.0463,
    0.0009,0.0012,0.0012,0.0008,0.0002,0.0011,0.0012,0.0008,0.0104,
    0.0027,0.0010,0.0007,0.0008,0.0011,0.0007,0.0010,0.0009,0.0004,0.0220,
    0.0036,0.0018,0.0011,0.0011,0.0014,0.0014,0.0015,0.0016,0.0008,0.0111,0.0442,
    0.0029,0.0061,0.0022,0.0020,0.0004,0.0028,0.0036,0.0020,0.0010,0.0012,0.0019,0.0190,
    0.0011,0.0006,0.0004,0.0003,0.0004,0.0007,0.0006,0.0005,0.0003,0.0025,0.0052,0.0007,0.0053,
    0.0014,0.0007,0.0006,0.0006,0.0005,0.0005,0.0006,0.0009,0.0007,0.0027,0.0052,0.0007,0.0010,0.0211,
    0.0021,0.0009,0.0007,0.0009,0.0003,0.0007,0.0012,0.0010,0.0004,0.0007,0.0012,0.0012,0.0003,0.0004,0.0221,
    0.0064,0.0020,0.0029,0.0024,0.0010,0.0017,0.0026,0.0034,0.0010,0.0015,0.0021,0.0026,0.0007,0.0010,0.0014,0.0167,
    0.0036,0.0015,0.0020,0.0016,0.0009,0.0012,0.0019,0.0019,0.0007,0.0024,0.0028,0.0020,0.0009,0.0011,0.0011,0.0048,0.0156,
    0.0003,0.0002,0.0001,0.0001,0.0001,0.0002,0.0002,0.0003,0.0001,0.0003,0.0006,0.0002,0.0002,0.0007,0.0001,0.0002,0.0002,0.0087,
    0.0011,0.0007,0.0006,0.0005,0.0003,0.0005,0.0006,0.0006,0.0016,0.0013,0.0020,0.0008,0.0005,0.0046,0.0003,0.0009,0.0008,0.0010,0.0148,
    0.0046,0.0013,0.0009,0.0010,0.0013,0.0010,0.0015,0.0014,0.0005,0.0123,0.0089,0.0015,0.0022,0.0022,0.0010,0.0021,0.0033,0.0004,0.0012,0.0246 };

}  // namnespace

namespace cs
{

BlosumMatrix::BlosumMatrix(Type matrix)
        : SubstitutionMatrix(AminoAcidAlphabet::instance()->size())
{
    switch (matrix) {
        case BLOSUM45:
            init(g_blosum45);
            break;
        case BLOSUM62:
            init(g_blosum62);
            break;
        case BLOSUM80:
            init(g_blosum80);
            break;
        default:
            throw MyException("Unsupported BLOSUM matrix!");
    }
}

void BlosumMatrix::init(const float* blosum_xx)
{
    // Read raw BLOSUM data vector
    int n = 0;
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b <= a; ++b, ++n)
            p_[a][b] = blosum_xx[n];

    // Add uppper right matrix part
    for (int a = 0; a < size_-1; ++a)
        for (int b = a+1; b < size_; ++b)
            p_[a][b] = p_[b][a];

    init_from_target_frequencies();
}

}  // cs
