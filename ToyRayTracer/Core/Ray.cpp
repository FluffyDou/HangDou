/******************************************************************/
/* Ray.cpp                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a couple methods for the Ray class.  Of       */
/*     particular interest may be the Ray::CheckHit() method,     */
/*     which simply checks if the current hit (at distance t)     */
/*     is closer than any previously hit surfaces.                */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "Ray.h"
//#include "DataTypes/MathDefs.h"

Ray::Ray()
{
	hitDist   = MAXFLOAT;
	hitObj    = 0;
    hitNormal = vec3(0.0, 0.0, 0.0);
    bounce    = 0;
}

Ray::Ray( const vec3 &orig, const vec3 &dir ): 
origin(orig), direction(dir)
{
    hitDist   = MAXFLOAT;
    hitObj    = 0;
    hitNormal = vec3(0.0, 0.0, 0.0);
    bounce    = 0;
}


Ray::Ray( const vec3 &orig, const vec3 &dir, const int bounce ): 
	origin(orig), direction(dir)
{
    hitDist   = MAXFLOAT;
    hitObj    = 0;
    hitNormal = vec3(0.0, 0.0, 0.0);
    this->bounce    = bounce;
}


int Ray::CheckHit( float t, Primitive *obj )
{
	// Do we miss?
	if (t >= hitDist || t <= EPSILON) return 0;

	// Nope.  A hit occurred
	hitDist = t;
	hitObj  = obj;

    return 1;
    
}



