/******************************************************************/
/* Sphere.h                                                       */
/* -----------------------                                        */
/*                                                                */
/* The file defines a sphere class that inherits from the basic   */
/*     object Primitive class.  This sphere renders a sphere      */
/*     of some radius centered about a given point.               */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef SPHERE_H
#define SPHERE_H

#include "Objects/Primitive.h"
#include "DataTypes/datatypes.h"

class Material;

class Sphere : public Primitive {
public:
	Sphere( Material *matl, const vec3 &center, float radius );
	virtual ~Sphere();
	virtual void Intersect( Ray &ray );
	virtual vec3 ComputeNormal( Ray &ray ) const;

private:
	vec3  center;
	float radius;
};

#endif

