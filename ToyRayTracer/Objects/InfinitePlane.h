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

#ifndef INFINITE_PLANE_H
#define INFINITE_PLANE_H

#include "Objects/Primitive.h"
#include "DataTypes/datatypes.h"

class Material;

class InfinitePlane : public Primitive {
public:
	InfinitePlane( Material *matl, const vec3 &norm, const vec3 &pointOnPlane );
	InfinitePlane( Material *matl, const vec4 &planeEq );
	virtual ~InfinitePlane();
	virtual void Intersect( Ray &ray );
	virtual vec3 ComputeNormal( Ray &ray ) const;

private:
	vec4 planeEq; // planeEq.X() is A, .Y() is B, .Z() is C, .W() is D
};

#endif

