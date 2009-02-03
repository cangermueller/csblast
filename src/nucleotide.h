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
    static inline Nucleotide& instance() { static Nucleotide inst; return inst; }

  protected:
    // Gets ctoi conversion array from derived class.
    virtual const char* get_itoc() const { return nucleotides_; }

  private:
    // IUPAC nucleotide code
    static const char nucleotides_[];

    Nucleotide();
    ~Nucleotide() {}
    // Disallow copy and assign
    Nucleotide(const Nucleotide& );
    Nucleotide& operator =(const Nucleotide& other);
};



const char Nucleotide::nucleotides_[] = "ACGT";

Nucleotide::Nucleotide()
        : Alphabet(4, 'N')
{ init(); }

}  // cs

#endif
