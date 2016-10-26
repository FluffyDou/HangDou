/*********************************************************************/
// LineSegment.h
// This the class holding line segments. It contains the two end
// points, mid point and the line direction (p0 -> p1)

// Hang 10/14/2014
/*********************************************************************/

#pragma once

#include "Eigen\Dense"

using Eigen::Vector3f;

class LineSegment
{
public:

    LineSegment(void);

    LineSegment(Vector3f p0, Vector3f p1);

    ~LineSegment(void);

    /*** Setters ***/
    void SetLineSeg(Vector3f p0, Vector3f p1);

    /*** Getters ***/
    // Return two end points
    const Vector3f& GetP0() const { return m_p0; }

    const Vector3f& GetP1() const { return m_p1; }

    // Return mid point
    const Vector3f& GetMidPoint() const { return m_midP; }

    // Return line direction
    const Vector3f& GetLineDir() const { return m_dir; }

    // Return line length
    const float& GetLen() const { return m_len; }

private:

    // Start point and end point
    Vector3f m_p0;
    Vector3f m_p1;
    // Line direction
    Vector3f m_dir;
    // Line segment mid point
    Vector3f m_midP;
    // Line length
    float    m_len;

};

