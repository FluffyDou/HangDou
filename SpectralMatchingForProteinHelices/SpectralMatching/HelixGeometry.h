/*********************************************************************/
// This the class holding helix geometry. Specifically, line segments.
// So far, it holds line segments and number of lines.
// The unit is Angstrom, the same as unit in PDB file (Mathematica use picometer).

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include <vector>
#include <array>
#include "Eigen\Dense"
#include "LineSegment.h"

using Eigen::Vector3f;

class HelixGeometry
{
public:
    typedef std::array<Vector3f,2> lineSeg;    

    // Create an empty object
    HelixGeometry(void);
    // Construct with allocating space for line segments
    HelixGeometry(int lineNum);

    ~HelixGeometry(void);

    /*** Setters ***/
    // Add a line segment into the list --- try to protect by m_lineNum. ToDo
    void AddLineSeg(int i, LineSegment lineSegment)  { m_lineSegList[i] = lineSegment; }
    // Add the residue abbreviation to sequence list
    void AddSeq(int i, const char abbrev)            { m_sequence[i].push_back(abbrev); }

    // Set model offset
    void SetModelOffset( Vector3f offset )           { m_modelOffset = offset; }

    // Compute the model center
    void ComputeModelCenter();

    // Transform the model to world space
    void TransformModel();

    /*** Getters ***/    
    // Get Line segment list -- Return the pointer to the lineSegs array
    const std::vector<LineSegment> GetLineSegList() const { return m_lineSegList; }

    // Get sequence list
    const std::vector< std::vector<char> > GetSeqList() const { return m_sequence; }

    // Return specific lineSeg --- try to protect by m_lineNum. ToDo
    const LineSegment& GetLineSegment(int i) const { return m_lineSegList[i]; }

    // Get model center and offset
    const Vector3f GetModelCenter() const { return m_modelCenter; }
    const Vector3f GetModelOffset() const { return m_modelOffset; }
    const int GetHelixNum() const         { return m_lineNum; }

    // Clear or release the contents
    void Clear();
    // Initialize --- Allocate memory
    void Reset(int lineNum);

private:
    
    // A list holding line segments
    // We use an vector of LineSegment to hold all the line segments
    std::vector<LineSegment> m_lineSegList;
    // A list holding sequence for this helix
    std::vector< std::vector<char> > m_sequence;

    // ATTENTION! Use these two only when dumping things onto the disk.
    // The model center
    Vector3f m_modelCenter;
    // Offset for visualize (in Mathematica)
    Vector3f m_modelOffset;

    // The number of line segments
    int m_lineNum;
};

