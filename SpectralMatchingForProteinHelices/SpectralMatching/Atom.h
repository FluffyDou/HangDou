/*********************************************************************/
// This the class holding all the data of an atom.
// So far, it holds coordinates only.

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include "Eigen\Dense"

using Eigen::Vector3f;

class Atom
{
public:
    // An empty atom
    Atom(void);
    // Construct an atom with coordinates
    Atom( Vector3f atomCoord );
    Atom(float x, float y, float z);

    ~Atom(void);

    /*** Setters ***/
    void SetCoord( Vector3f coord )            { m_coord = coord; }
    void SetCoord( float x, float y, float z ) { m_coord = Vector3f(x, y, z); }

    /*** Getters ***/
    const Vector3f& GetCoord() const { return m_coord; }

private:

    // 3D Coordinate of the atom
    Vector3f m_coord;
};

