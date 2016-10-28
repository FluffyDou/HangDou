/******************************************************************/
/* Cylinder.h                                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines a cylinder class that inherits from the basic */
/*     object Primitive class.  This class renders a cylinder     */
/*     based on a radius, height, axis, and starting point.       */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef CYLINDER_H
#define CYLINDER_H

#include "Objects/Primitive.h"
#include "DataTypes/datatypes.h"

class Material;

class Cylinder : public Primitive {
public:
	Cylinder( Material *matl, const vec3 &center, const vec3 &axis, float radius, float height );
	virtual ~Cylinder();
	virtual void Intersect( Ray &ray );
	virtual vec3 ComputeNormal( Ray &ray ) const;

private:
	// Internal memory to store cylinder position & orientation
	vec3 center, axis;
	float radius, height;

	// Internal memory for reusable class computations
	vec3 uVec, vVec;
};

#endif

