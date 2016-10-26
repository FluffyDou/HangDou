/*********************************************************************/
// AlphaHelix.h
// This is a class holding alpha helix. It holds the start and end index 
// of the residues constructing this helix.
// ATTENTION! Here we hold offset indices. Use them to directly access residues.

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include <vector>

class AlphaHelix
{
public:

    AlphaHelix(void);

    AlphaHelix( std::vector< int > indices );

    ~AlphaHelix(void);

    /******* Setters ******/
    // Add residue indices (indices should be a vec2: start and end index)
    void AddResidueIndices( std::vector< int > indices )  { m_residueIndices = indices; }

    // Set residue indices --- add offset
    void SetResidueIndices( std::vector< int > indices )  { m_residueIndices = indices; }

    /******* Getters ******/
    // Return the start and end residue indices of this helix
    const std::vector< int >& GetResidueIndices() const { return m_residueIndices; }

private:
    // The list or vector holding all the residues for this helix
    // Todo: change to  std::array<int,2>
    std::vector< int > m_residueIndices;    
};

