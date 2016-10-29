/******************************************************************/
/* igluHalton.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* The file defines a quasi random number generator using the     */
/*    Halton sequence.                                            */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IGLUHALTON_H
#define IGLUHALTON_H

#include <time.h>

namespace iglu {

class IGLUHalton
{
public:
	IGLUHalton( unsigned int baseX=2, unsigned int baseY=3, unsigned int baseZ=5 );
	~IGLUHalton() {}

	// Really, only the 
	double         dHalton( unsigned int idx );          /* in [0..1]     */
	float          fHalton( unsigned int idx );          /* in [0..1]     */
	
	/* Give a random vector on the unit sphere, +z unit hemisphere, or arbitrary unit hemisphere */
#if 0
	vec3 RandomSphereVector( void );                    
	vec3 RandomHemisphereVector( void );               
	vec3 RandomHemisphereVector( const vec3 &inDir ); 
	vec3 RandomCosineWeightedHemisphereVector( const vec3 &inDir ); 
#endif
};

}


#endif

