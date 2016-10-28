/******************************************************************/
/* Sphere.cpp                                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines a sphere class that inherits from the basic   */
/*     object Primitive class.  This sphere renders a sphere      */
/*     of some radius centered about a given point.               */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "Core/Ray.h"
#include "Core/Scene.h"
#include "Objects/Sphere.h"
#include "DataTypes/MathDefs.h"
#include <math.h>

Sphere::Sphere( Material *matl, const vec3 &center, float radius ) :
	Primitive(matl)
{
	this->center = center;
    this->radius = radius;
}

Sphere::~Sphere()
{
	// Fill in sphere destructor (if needed)
}


void Sphere::Intersect( Ray &ray ) 
{
    vec3  EO = center - ray.origin;
    float v  = EO.Dot( ray.direction );
    float d  = radius*radius - ( EO.Dot(EO) - v*v);

    if(d >= EPSILON)
    {
        if(ray.CheckHit( v - sqrtf(d), this ) > 0)
            ray.hitNormal = this->ComputeNormal( ray );
    }
}


vec3 Sphere::ComputeNormal( Ray &ray ) const
{
	// Fill in sphere normal computations
    return ( ray.GetHitPoint() - center ).vNormalize();
}

