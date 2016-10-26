/*********************************************************************/
// AminoAcidResidue.h
// This the class holding all the data of a amino acid residue.
// So far, it holds atoms, amino abbreviation only.

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include <vector>
#include "Atom.h"

class AminoAcidResidue
{
public:
    typedef std::vector<Atom> atomList;

    // Empty constructor
    AminoAcidResidue(void);
    // Construct a residue with amino acid abbreviation
    AminoAcidResidue(char aminoAbbrev);

    ~AminoAcidResidue(void);

    /******* Setters ******/
    void AddOneAtom( Atom atom ) { m_atoms.push_back(atom); }

    /******* Getters ******/
    // Get all the atoms for this residue
    const atomList& GetAtoms()    const { return m_atoms; }
    // Get c-alpha atom for this residue
    const Atom& GetCAlpha()       const { return m_atoms[C_ALPHA_LOC]; }
    // Get amino abbreviation character
    const char GetAminoAbbrev()  const { return m_aminoAbbrev; }

private:

    // Atoms in this residue
    atomList m_atoms;

    // Here we assume the second atom in the residue is the c-alpha carbon
    static const int C_ALPHA_LOC = 1;

    // Amino abbreviation
    char m_aminoAbbrev;
};

