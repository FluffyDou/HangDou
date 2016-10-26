/*********************************************************************/
// ProteinChain.h
// This the class holding all the data of a protein chain.
// Specifically, it holds (picked) information of a picked chain
// in a standard PDB file (.pdb), like atom coordinates, sequence,
// secondary structure sequence and so on

// Chain:
//      Sequence:
//              char list of Amino Abbreviation
//      SSE:
//              Helix: a map holding residue index
//      Residues:
//              char of Amino Abbreviation
//              Atoms:
//                      Atom:
//                              3D point coordinate.

// Hang 10/06/2014
/*********************************************************************/

#pragma once

#include <vector>
#include <map>
#include "AminoAcidResidue.h"
#include "SecondaryStructure.h"
#include "AlphaHelix.h"

class ProteinChain
{
public:

    typedef AminoAcidResidue Residue;
    typedef std::vector<Residue> residueList;
    typedef std::map<int, Residue> residuePool;
    typedef std::vector<char> charList;

	// Create an empty chain
	ProteinChain(void);

    /******* Setters ******/
    void SetChainID(char id)         { m_chainID = id; }
    //void SetResidueNum(int num)    { m_residueNum = num; }

    // Add residue abbreviation into the list
    void AddSeqChar(char seqChar)    { m_sequence.push_back(seqChar); }
    // Add a residue into the list
    void AddResidue(int index, Residue residue) { m_residues[index] = residue; }

    // Add an atom to one residue
    void AddAtomToResidue(int residueIndex, Atom atom)  
    { 
        // In PDB file, residueIndex may not start from 1. 
        // Should protect on "residueIndex". ToDo.
        m_residues[residueIndex].AddOneAtom(atom);
    }

    // Add Helix to the helix list
    void AddHelixToSSE( AlphaHelix helix ) { m_sse.AddHelix(helix); }

    /******* Getters ******/
    const char GetChainID()   const  { return m_chainID; }    

    const charList&     GetSequence() const { return m_sequence; }
    const residuePool&  GetAllResidues() const { return m_residues; }

    // Check if memory has been allocated for stuffs in this chain
	//bool IsInitialized()  { return m_alreadyInit; }

    // Mainly to add offset to alpha helix (residue index)
    SecondaryStructure& GetSSEToUpdate() { return m_sse; }

    // Return the secondary structure elements
    const SecondaryStructure& GetSSE() const { return m_sse; }

	~ProteinChain(void);

private:

    // Valid ID is from 'A' to 'Y'.
    char  m_chainID;
    int   m_residueNum;
	//bool  m_alreadyInit;

    // Hold the sequence character for the chain
    charList m_sequence;

    // A map holding all the residues for the chain    
    residuePool m_residues;

    // Hold all the secondary structure elements
    SecondaryStructure m_sse;
    
};
