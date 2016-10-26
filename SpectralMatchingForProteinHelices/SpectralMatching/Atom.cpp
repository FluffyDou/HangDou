/*********************************************************************/
// This the class holding all the data of an atom.
// So far, it holds coordinates only.

// Hang 10/08/2014
/*********************************************************************/

#include "Atom.h"

// Construct an empty atom
Atom::Atom(void)
{
    m_coord = Vector3f(0, 0, 0);
}

// Construct an atom with coordinates
Atom::Atom(float x, float y, float z)
{
    m_coord = Vector3f(x, y, z);
}

Atom::Atom( Vector3f atomCoord )
{
    m_coord = atomCoord;
}

Atom::~Atom(void)
{
}
