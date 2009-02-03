#ifndef CS_AMINO_ACID_H
#define CS_AMINO_ACID_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about amino acids

#include "alphabet.h"

namespace cs
{

class AminoAcid : public Alphabet
{
  public:
    static inline AminoAcid& instance() { static AminoAcid inst; return inst; }

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return amino_acids_; }

  private:
    // IUPAC amino acid code
    static const char amino_acids_[];

    AminoAcid();
    ~AminoAcid() {}
    // Disallow copy and assign
    AminoAcid(const AminoAcid& );
    AminoAcid& operator =(const AminoAcid& other);
};


const char AminoAcid::amino_acids_[] = "ARNDCQEGHILKMFPSTWYV";

AminoAcid::AminoAcid()
        : Alphabet(20, 'X')
{
    init();
    // Conversion of non-standard amino acids to integers
    set_ctoi('O', any());
    set_ctoi('J', any());
    set_ctoi('U', ctoi('C'));
    set_ctoi('B', ctoi('D'));
    set_ctoi('Z', ctoi('E'));
}

}  // cs

#endif
