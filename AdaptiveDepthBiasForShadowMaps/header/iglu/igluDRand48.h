/******************************************************************/
/* igluRandom.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* The file defines a random number generator class.  This may    */
/*   not be the best random number generator around, but it works */
/*   reasonably well, and gives a portable random class for use   */
/*   on multiple platforms (which may or may not have rand() or   */
/*   drand48(), or they may be broken).                           */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef IGLUDRAND48_H
#define IGLUDRAND48_H

#include <time.h>
#include <stdint.h>

namespace iglu {

class vec3;

class IGLUDRand48
{
private:
	ushort dRandSeeds[3];
	ushort dRandMultipler[3]; // a 48-bit multiplier
	ulong  dRandAdder;

	void DoRandomStep( void );

public:
	IGLUDRand48( uint32_t seed = (ulong)time(0) );
	IGLUDRand48( uint64_t seed );
	~IGLUDRand48() {}

	double         dRandom( void );          /* in [0..1]      */
	float          fRandom( void );          /* in [0..1]      */
	unsigned short sRandom( void );          /* in [0..2^16-1] */
	unsigned int   iRandom( void );          /* in [0..2^32-1] */
	bool           bRandom( void );          /* true or false  */
	unsigned char  cRandom( void );          /* in [0..255]    */
	
	/* Give a random vector on the unit sphere, +z unit hemisphere, or arbitrary unit hemisphere */
	vec3 RandomSphereVector( void );                    
	vec3 RandomHemisphereVector( void );               
	vec3 RandomHemisphereVector( const vec3 &inDir ); 
	vec3 RandomCosineWeightedHemisphereVector( const vec3 &inDir ); 

	// A pointer to a IGLURandom could have type IGLURandom::Ptr
	typedef IGLUDRand48 *Ptr;
};

}


#endif

