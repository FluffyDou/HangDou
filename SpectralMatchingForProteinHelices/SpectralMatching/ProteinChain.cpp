/*********************************************************************/
// ProteinChain.cpp
// This the class holding all the data of a protein chain.
// Specifically, it holds (picked) information of a picked chain
// in a standard PDB file (.pdb), like atom coordinates, sequence,
//  secondary structure sequence and so on

// Hang 10/06/2014
/*********************************************************************/

#include "ProteinChain.h"


ProteinChain::ProteinChain(void)
{
	//m_alreadyInit = false;	
    //m_seqOffset   = 0;
	m_chainID     = 0;
	m_residueNum  = 0;
}


//// So far Deprecated since we use vectors
//void ProteinChain::InitContentIfNotYet()
//{ 
//	// Check if initialized before
//	if( m_alreadyInit )
//		return;	
//
//    // We don't allocate memory for sequence vector
//    // because of sequence reading convenience in PDB
//
//
//	m_alreadyInit = true;
//}

ProteinChain::~ProteinChain(void)
{
    //if( 0 != m_sequence )
	//    delete[] m_sequence;
}

