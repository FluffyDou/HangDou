/******************************************************************/
/* InfinitePlane.h                                                */
/* -----------------------                                        */
/*                                                                */
/* The file defines an infinite plane class that inherits from    */
/*     the basic object Primitive class.  This class renders an   */
/*     infinite plane solving the equation Ax+By+Cz+D=0.  Here,   */
/*     (A,B,C) represents the unit surface normal from the plane  */
/*     and D is the distance from the plane to the origin along   */
/*     the vector (A,B,C).   One can also compute D based upon    */
/*     the surface normal (A,B,C) and a point P=(x,y,z) on the    */
/*     plane.                                                     */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/


#include "Core/Ray.h"
#include "Core/Scene.h"
#include "Objects/InfinitePlane.h"
#include "DataTypes/MathDefs.h"
#include <math.h>

InfinitePlane::InfinitePlane( Material *matl, const vec3 &norm, const vec3 &pointOnPlane ) :
	Primitive(matl)
{
	// Fill in constructor #1
}

InfinitePlane::InfinitePlane( Material *matl, const vec4 &planeEq ) :
	Primitive(matl)
{
	// Fill in constructor #2
}

InfinitePlane::~InfinitePlane()
{
	// Fill in plane destructor (if needed)
}


void InfinitePlane::Intersect( Ray &ray ) 
{
	// Fill in plane intersector
}


vec3 InfinitePlane::ComputeNormal( Ray &ray ) const
{
	// Fill in plane normal compuations
    return vec3(0.0, 0.0, 0.0);
}

