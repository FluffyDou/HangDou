/******************************************************************/
/* Object.h                                                       */
/* -----------------------                                        */
/*                                                                */
/* The file defines a base object class that is inherited by more */
/*     complex object types.                                      */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef OBJECT_H
#define OBJECT_H

#include <stdio.h>

class Ray;
class Material;
class Vector;
class Scene;

class Object {
public:
	Object() {}
	virtual ~Object() {}
	virtual void Intersect( Ray &ray ) = 0;
};


#endif

