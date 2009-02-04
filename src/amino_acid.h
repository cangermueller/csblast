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
    static const AminoAcid& instance();

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return amino_acids_; }

  private:
    // IUPAC amino acid code
    static const char amino_acids_[];

    AminoAcid();
    virtual ~AminoAcid() {}

    // Disallow copy and assign
    AminoAcid(const AminoAcid& );
    AminoAcid& operator =(const AminoAcid& other);
};

}  // cs

#endif
