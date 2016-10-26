/*********************************************************************/
// LineSegment.cpp
// This the class holding line segments. It contains the two end
// points, mid point and the line direction (p0 -> p1)

// Hang 10/14/2014
/*********************************************************************/

#include "LineSegment.h"

LineSegment::LineSegment(void)
{
    Vector3f temp(0,0,0);
    m_p0   = temp;
    m_p1   = temp;
    // Line direction
    m_dir  = temp;
    // Line segment mid point
    m_midP = temp;
}


LineSegment::LineSegment(Vector3f p0, Vector3f p1)
{
    SetLineSeg(p0, p1);
}


LineSegment::~LineSegment(void)
{
}

void LineSegment::SetLineSeg(Vector3f p0, Vector3f p1)
{
    // Record two end points
    m_p0   = p0;
    m_p1   = p1;
    // Record line direction and mid point
    m_dir  = (p1 - p0).normalized();
    m_midP = 0.5*(p0 + p1);
    // Record the line segment length
    m_len  = (p1 - p0).norm();
}