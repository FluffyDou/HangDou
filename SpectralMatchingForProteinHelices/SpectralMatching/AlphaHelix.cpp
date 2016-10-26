/*********************************************************************/
// AlphaHelix.h
// This is a class holding alpha helices.
// So far, it holds a list of residue indices of start and end residue 
// on the protein chain

// Hang 10/08/2014
/*********************************************************************/

#include "AlphaHelix.h"

AlphaHelix::AlphaHelix(void)
{
}

// This constructor should always take vec2
 AlphaHelix::AlphaHelix(std::vector< int > indices)
 {
     m_residueIndices = indices;
 }

AlphaHelix::~AlphaHelix(void)
{
}
