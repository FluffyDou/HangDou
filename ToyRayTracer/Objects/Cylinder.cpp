/******************************************************************/
/* Cylinder.cpp                                                   */
/* -----------------------                                        */
/*                                                                */
/* The file defines a cylinder class that inherits from the basic */
/*     object Primitive class.  This class renders a cylinder     */
/*     based on a radius, height, axis, and starting point.       */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "Core/Ray.h"
#include "Core/Scene.h"
#include "Objects/Cylinder.h"
#include "DataTypes/MathDefs.h"
#include <math.h>

Cylinder::Cylinder( Material *matl, const vec3 &center, const vec3 &axis, float radius, float height ) :
	Primitive(matl), // The cylinder has an internal, public Primitive class.  Initialize that.
	center(center),  // Copy the input to the constructor into the internal object storage
	axis(axis),      //   same for the input axis
 	radius(radius),  //   same for the input radius
	height(height)   //   same for the input height
{
	// Now, compute two vector (uVec and vVec) that are perpendicular to the cylinder
	//    axis as well as being mutually perpendicular (i.e., compute an ortho normal basis)
	vec3 vec = axis.Cross(vec3::XAxis());    // Cross cyl axis by x-axis to get a vector perp to both
	if (vec.LengthSqr() < 0.01)              // Is it too short?  =>  Numerical errors
		vec = axis.Cross(vec3::YAxis());     //    If too short, compute a new vec perp to cyl axis & y-axis
	uVec = axis.Cross( vec ).vNormalize();   // Do some cross products & normalization to find two vectors
	vVec = axis.Cross( uVec ).vNormalize();  //    perpendicular to each other & the axis.
}

Cylinder::~Cylinder()
{
	// The cylinder does not allocate memory, so it need not deallocate any
}

// This does intersection for a cylinder.  I'm not going to comment this well, since you can
//    run the math yourself to figure out what is going on.  The key is to understand that
//    some math happens, you compute zero, one, or more candidate hits, then you call 
//    ray.CheckHit() which guarantees the ray remembers the geometry closest to its origin.
void Cylinder::Intersect( Ray &ray ) 
{
	vec3 tmp  = ray.origin-center;
	vec3 orig = vec3( tmp.Dot( uVec ), tmp.Dot( vVec ), tmp.Dot( axis ) );
	vec3 dir  = vec3( ray.direction.Dot( uVec ), ray.direction.Dot( vVec ), ray.direction.Dot( axis ) );

	float a = dir.X()*dir.X()+dir.Y()*dir.Y();
	if (a <= 0) return;  // If the ray is parallel to the cylinder we miss.

	float b = 2*(orig.X()*dir.X()+orig.Y()*dir.Y());
	float c = orig.X()*orig.X()+orig.Y()*orig.Y()-radius*radius;
	float d = b*b-4*a*c;
	if (d < 0) return;  // We miss.  (Radical in quadratic eq is < 0)
	
	d = sqrtf(d);
	a *= 2;
	float t1 = (-b+d)/a;
	float t2 = (-b-d)/a;
	if (t1>t2) // if ordered wrong (so t1 is further than t2), swap.
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
	}
	float z1 = orig.Z()+t1*dir.Z();
	float z2 = orig.Z()+t2*dir.Z();

	// We have two candidate hitpoints.  One a distance t1 from the ray origin,
	//    one a distance t2.
	// In order for these to hit, they have to be on the cylinder surface at an
	//    appropriate height (i.e., if they're too far from the center of the 
	//    cylinder, they don't hit)
	if (t1 > 0 && z1 >= -0.5*height && z1 <= 0.5*height)
    {
        if(ray.CheckHit( t1, this ) != 0)
            ray.hitNormal = this->ComputeNormal( ray );
		
    }
	else if (t2 > 0 && z2 >= -0.5*height && z2 <= 0.5*height)
    {
        if( ray.CheckHit( t2, this ) != 0 )
            ray.hitNormal = this->ComputeNormal( ray );
    }
}

// The surface normal of a cylinder points outwards radially from the axis.
//    1) So first find the hit point.
//    2) Find the vector from the "center" of the cylinder radially outwards
//       through the point on the surface.
//    3) This vector radiates outwards, but also has traverses down the cylinder
//       axis.  Subtract out the component that points in the direction of the
//       axis, giving just the vector pointing radially outwards.
//    4) Normalize the result and return it.
vec3 Cylinder::ComputeNormal( Ray &ray ) const
{
	vec3 vec = ray.GetHitPoint()-center;
	float dotPrd = vec.Dot(axis);
	vec = vec + (-dotPrd)*axis;
	return vec.vNormalize();
}

