/*********************************************************************/
// PDBObject.h
// This is a class holding all the data from standard PDB file.

// A PDB object has a parser to read in standard PDB file (.pdb).
// A PDB object holds a list of Chain Object, ProteinChain.

// A ProteinChain has information as that in Mathematica, like atom 
// coordinates, sequence, secondary structure sequence and so on

// Hang 10/06/2014
/*********************************************************************/

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>

#include "ProteinChain.h"

class PDBObject
{
public:

    typedef AminoAcidResidue Residue;
    typedef std::vector<ProteinChain> proteinChainList;

    // Get the abbreviation of an amino acid
    static char GetAminoAbbrev(const char * aminoAcid);

	// Create an empty object
	PDBObject(void);
	// Read in all the chains from the input .pdb file.
	PDBObject::PDBObject(char* fileName);
	// Read in specific chain information from the input .pdb file.
	PDBObject::PDBObject(char* fileName, int chainID);
    
    /******* Setters ******/

    /******* Getters ******/
    // Return all the chains
    //const proteinChainList& GetChains() const { return m_chains; }
    const std::map<char, ProteinChain> & GetChains() const { return m_proteinChains; }

    // Return a specific chain --- protect by index using a map. ToDo.
    //const ProteinChain& GetChain(char i) const { return m_chains[i - m_chainIdOffset]; }
    const ProteinChain GetChain(char i) const
    { 
        auto chain = m_proteinChains.find(i);
        if(chain != m_proteinChains.end()) 
            return chain->second;

        return ProteinChain();
    }
    // Return chain offset
    //const char GetChainIdOffset() const { return m_chainIdOffset; }

    /****** Others ******/
    // Read in a PDB file, only supposed to be called by the constructor
    void ReadInPDBFile(const char* fileName);

    // Clear the content
    //void Clear() { m_chains.clear(); m_proteinChains.clear(); }
    void Clear() { m_proteinChains.clear(); }

    typedef PDBObject *Ptr;

	~PDBObject(void);

private:
     // Retrieve information from the read in line from PDB file
     int         readInt(std::string line, int start, int end);
     float       readFloat(std::string line, int start, int end);
     std::string readString(std::string line, int start, int end);
     // Trim the space in the given string
     std::string trimString(std::string const& str);

private:

    //proteinChainList m_chains;
    //char m_chainIdOffset;

    std::map<char, ProteinChain> m_proteinChains;
};

