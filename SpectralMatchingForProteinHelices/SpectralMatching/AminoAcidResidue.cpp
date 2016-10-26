/*********************************************************************/
// AminoAcidResidue.cpp
// This the class holding all the data of a amino acid residue.
// So far, it holds atoms only.

// Hang 10/08/2014
/*********************************************************************/

#include "AminoAcidResidue.h"


AminoAcidResidue::AminoAcidResidue(void)
{
    m_aminoAbbrev = 0;
}

// Construct a residue with amino acid abbreviation
AminoAcidResidue::AminoAcidResidue(char aminoAbbrev)
{
    m_aminoAbbrev = aminoAbbrev;
}

AminoAcidResidue::~AminoAcidResidue(void)
{
}
