/******************************************************************/
/* Ray.h                                                          */
/* -----------------------                                        */
/*                                                                */
/* The file defines a ray class that stores relevent information  */
/*     to a ray.  This includes its origin, its direction, and    */
/*     the current ray depth (how many times it has bounced).     */
/*     It also stores hit information (which is the closest       */
/*     object hit so far).                                        */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef RAY_H
#define RAY_H 1

#include "DataTypes/vec3.h"
#include "DataTypes/MathDefs.h"

class Object;
class Primitive;

// This "class" is really a glorified structure.  

class Ray
{
public:  // Data

  // basic information about a ray
  vec3 origin;
  vec3 direction;

  // Information about the hit point
  float hitDist;
  Primitive *hitObj;

  // record the hit primitive
  vec3 hitNormal;

  // record the bounce number
  int bounce;

public:  // Methods
	Ray();

	// Set up the ray with a give origin and direction (and optionally a ray depth)
    Ray( const vec3 &orig, const vec3 &dir);
	Ray( const vec3 &orig, const vec3 &dir, const int bounce );

	// Call this method when you hit an object obj at a distance t from the ray origin.
	//    It tests to see if t is the closest object hit yet, and if so remembers it.
    // When hit something new return 1, while when miss, return 0
	int CheckHit( float t, Primitive *obj );
	
	// Gets the hit point and/or object.  If no point was hit, the HitObject is NULL.
	inline vec3 GetHitPoint( void ) const { return origin + hitDist*direction; }

	// Checks if there has been an intersection along this ray.
	inline bool WasIntersection( void ) const { return hitObj != 0; }
};


#endif

