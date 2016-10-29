/******************************************************************/
/* igluSampler2D.h                                                */
/* -----------------------                                        */
/*                                                                */
/* The file defines a base class from which all 2D samplers are   */
/*     derived.  Basically, base classes have to define how to    */
/*     get a sample in the range x,y in [0..1), and this base     */
/*     class defines methods that can use these samples to get    */
/*     more interesting behaviors (e.g., samples on hemispheres)  */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUSAMPLER2D_H
#define IGLUSAMPLER2D_H

namespace iglu {

class vec2;
class vec3;

class IGLUSampler2D
{
public:
	IGLUSampler2D( void ) {}
	~IGLUSampler2D() {}

	///////////////////////////////////////////////////////////////////////////
	// This is the function all derived classes must implement.
	//    It needs to give a sample for x,y in [0..1)
	virtual vec2 Sample( unsigned int sampIdx = 0 ) = 0;


	///////////////////////////////////////////////////////////////////////////
	// All functions below use Sample().  Because of that, they need not be 
	//    redefined in derived classes.  Only Sample() needs definition!
	//    All require your current sample number.  For pure random sampling,
	//    you may be able to ignore the sampIdx and always use 0.
	
	///////////////////////////////////////////////////////////////////////////
	// These sample 2D surfaces (so are clearly 2D)

	// Sample the unit disk.  Return values will satisfy (x^2 + y^2 <= 1)
	vec2 SampleDisk( unsigned int sampIdx = 0 );

	// Sample the triangle with specified vertices uniformly
	vec3 SampleTriangle( const vec3 &v0, const vec3 &v1, const vec3 &v2, unsigned int sampIdx = 0 );

	///////////////////////////////////////////////////////////////////////////
	// These give 3D values, but they are inherently 2D, since they are surfaces.

	// Sample the unit sphere uniformly (x^2 + y^2 + z^2 = 1 in return value)
	vec3 SampleSphere( unsigned int sampIdx = 0 );

    // Sample the unit +z hemisphere uniformly (x^2 + y^2 + z^2 = 1, z >= 0)
	vec3 SamplePositiveZHemisphere( unsigned int sampIdx = 0 );             

	// Sample the unit +z hemisphere with a cos distribution (more samples near x,y=0)
	vec3 SampleCosineWeightedPositiveZHemisphere( unsigned int sampIdx = 0 );        

	// Sample an arbitrary hemisphere oriented in a particular direction
	vec3 SampleHemisphere( const vec3 &inDir, unsigned int sampIdx = 0 ); 
	vec3 SampleCosineWeightedHemisphere( const vec3 &inDir, unsigned int sampIdx = 0 ); 

	// A pointer to a IGLUSampler2D could have type IGLUSampler2D::Ptr
	typedef IGLUSampler2D *Ptr;

};


// End iglu namespace
}


#endif

