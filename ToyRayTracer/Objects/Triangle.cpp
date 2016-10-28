/******************************************************************/
/* Triangle.cpp                                                   */
/* -----------------------                                        */
/*                                                                */
/* The file defines a triangle class that inherits from the basic */
/*     object Primitive class.  This sphere renders a triangle    */
/*     based on three vertex locations.                           */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/


#include "Core/Ray.h"
#include "Core/Scene.h"
#include "Objects/Triangle.h"
#include "DataTypes/MathDefs.h"
#include <math.h>

Triangle::Triangle( Material *matl, const vec3 &v0, const vec3 &v1, const vec3 &v2 ) :
	Primitive(matl), v0(v0), v1(v1), v2(v2)
{
	// Fill in constructor
}

Triangle::~Triangle()
{
	// Fill in triangle destructor (if needed)
}


void Triangle::Intersect( Ray &ray ) 
{
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v0 - v2;
    vec3 norm  = edge1.Cross(edge2);

    // Find intersection on the plane
    float r   = 1.0f / norm.Dot( ray.direction );
    vec3  e2  = v0 - ray.origin;
    float va  = norm.Dot( e2 );
    float t   = r*va;

    // Find the barycentric coordinates
    vec3  i  = e2.Cross( ray.direction );
    float u  = r * i.Dot( edge2 );
    float v  = r * i.Dot( edge1 );

    // If u & v are valid, check the hit point
    if( (u >= EPSILON) && (v >= EPSILON) && (u+v <= 1.0f)) 
    {
        if( ray.CheckHit( t, this ) != 0 )
            ray.hitNormal = this->ComputeNormal( ray );
    }

}


vec3 Triangle::ComputeNormal( Ray &ray ) const
{
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v0 - v2;
    return ( edge1.Cross(edge2) ).vNormalize();
}

