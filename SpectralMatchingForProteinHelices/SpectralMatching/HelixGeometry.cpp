/*********************************************************************/
// This the class holding helix geometry. Specifically, line segments.
// So far, it holds line segments and number of lines.
// The unit is Angstrom, the same as unit in PDB file (Mathematica use picometer).

// Hang 10/08/2014
/*********************************************************************/

#include "HelixGeometry.h"

// Create an empty object
HelixGeometry::HelixGeometry(void)
{    
    Clear();
}

// Initialize or reset --- Allocate memory
void HelixGeometry::Reset(int lineNum)
{
    // Clear old stuff
    Clear();
    // Set up new stuff
    m_lineNum  = lineNum;

    m_lineSegList = std::vector<LineSegment>(lineNum);
    m_sequence    = std::vector< std::vector<char> >(lineNum);
}

// Construct with allocating space for line segments
HelixGeometry::HelixGeometry(int lineNum)
{    
    Reset(lineNum);
}

HelixGeometry::~HelixGeometry(void)
{
}

void HelixGeometry::Clear()
{
    m_sequence.clear();
    m_lineSegList.clear();

    m_lineNum     = 0;
    m_modelCenter = Vector3f(0,0,0);
    m_modelOffset = Vector3f(0,0,0);    
}

// Compute the model center
void HelixGeometry::ComputeModelCenter()
{
    for( auto it = m_lineSegList.begin(); it != m_lineSegList.end(); ++it )
    {
        m_modelCenter += it->GetP0() + it->GetP1();
    }

    m_modelCenter /= (float)(2*m_lineNum);
}

// Transform the model to world space
void HelixGeometry::TransformModel()
{
    for( auto it = m_lineSegList.begin(); it != m_lineSegList.end(); ++it )
    {
        it->SetLineSeg( it->GetP0() - m_modelCenter + m_modelOffset, it->GetP1() - m_modelCenter + m_modelOffset );
    }
}