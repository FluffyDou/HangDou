/*********************************************************************/
// SecondaryStructure.h
// This is a class holding secondary structure of a protein chain.
// So far, it holds alpha helix list.

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include <vector>
#include "AlphaHelix.h"

class SecondaryStructure
{
public:

    //typedef std::vector<AlphaHelix> helixList;

    SecondaryStructure(void);

    ~SecondaryStructure(void);

    /******* Setters ******/
    // Add a helix to the helix list
    void AddHelix(AlphaHelix helix) { m_helices.push_back(helix); }

    /******* Getters ******/
    // Return the helix list
    const std::vector<AlphaHelix>& GetHelices() const { return m_helices; }
    // Return the length of helix list or helix number
    const int GetHeliceNum() const { return m_helices.size(); }

private:

    // A list of alpha helix
   std::vector<AlphaHelix> m_helices;   
};

