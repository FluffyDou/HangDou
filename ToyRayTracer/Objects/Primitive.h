/******************************************************************/
/* Primitive.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* The file defines a primitive (aka drawable) geometry class     */
/*     that is inherited by more complex objects (e.g., spheres,  */
/*     triangles, planes, etc).  Note this is different than a    */
/*     simple Object class, as it stores material properties,     */
/*     and allows the determination of a surface normal, both of  */
/*     which are not necessary for container objects, like Group. */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "Objects/Object.h"

class Ray;
class Material;
class Scene;

class Primitive : public Object {
public:
	Primitive(Material *matl=0) : matl(matl) {}
	virtual ~Primitive() {}
	virtual void Intersect( Ray &ray ) = 0;
	virtual vec3 ComputeNormal( Ray &ray ) const = 0;
	inline Material *GetMaterial( void ) { return matl; }
private:
	Material *matl;
};


#endif

