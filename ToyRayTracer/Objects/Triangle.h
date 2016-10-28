/******************************************************************/
/* Triangle.h                                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines a triangle class that inherits from the basic */
/*     object Primitive class.  This sphere renders a triangle    */
/*     based on three vertex locations.                           */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Objects/Primitive.h"
#include "DataTypes/datatypes.h"

class Material;

class Triangle : public Primitive {
public:
	Triangle( Material *matl, const vec3 &v0, const vec3 &v1, const vec3 &v2 );
	virtual ~Triangle();
	virtual void Intersect( Ray &ray );
	virtual vec3 ComputeNormal( Ray &ray ) const;

private:
	vec3 v0, v1, v2;
};

#endif

