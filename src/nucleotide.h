#ifndef CS_NUCLEOTIDE_H
#define CS_NUCLEOTIDE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Singleton class that encapsulates meta information about nucleic acids

#include "alphabet.h"

namespace cs
{

class Nucleotide : public Alphabet
{
  public:
    static Nucleotide& instance();

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return nucleotides_; }

  private:
    // IUPAC nucleotide code
    static const char nucleotides_[];

    Nucleotide();
    virtual ~Nucleotide() {}

    // Disallow copy and assign
    Nucleotide(const Nucleotide& );
    Nucleotide& operator =(const Nucleotide& other);
};

}  // cs

#endif
