/******************************************************************/
/* Random.h                                                       */
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
#include <math.h>
#include "Utils/Random.h"

Random::Random( unsigned long seed ) :
   maxshort( 65536L ), multiplier( 1194211693L ), adder( 12345L )
{
	randSeed = ( seed ? seed : (unsigned long)time(0) );
}

double Random::dRandom( void )
{
	randSeed = multiplier * randSeed + adder;
	return ((randSeed >> 16) % maxshort) / (double)maxshort;
}

float  Random::fRandom( void )
{
	randSeed = multiplier * randSeed + adder;
	return ((randSeed >> 16) % maxshort) / (float)maxshort;
}

unsigned short Random::sRandom ( void )
{
	randSeed = multiplier * randSeed + adder;
	return (unsigned short)((randSeed >> 16) % maxshort);
}

unsigned char Random::cRandom ( void )
{
	randSeed = multiplier * randSeed + adder;
	return (unsigned char)((randSeed >> 16) % 256L );
}

bool Random::bRandom ( void )
{
	randSeed = multiplier * randSeed + adder;
	return ((randSeed >> 16) & 0x00000001) > 0;  // check if the last bit is on or off
}

vec3 Random::RandomHemisphereVector( void )
{
	float f1 = fRandom(), f2 = fRandom();
	return vec3( sqrt(f1)*cos(2.0f*M_PI*f2), 
				   sqrt(f1)*sin(2.0f*M_PI*f2), 
				     sqrt(1.0f - f1) );
}

vec3 Random::RandomSphereVector( void )
{
	float f1 = 2.0f * fRandom() - 1.0f;
	float f2 = 2.0f * (float)M_PI * fRandom();
	float r = sqrt( 1.0f - f1*f1 );
	return vec3( r*cos( f2 ), r*sin( f2 ), f1 );
}

vec3 Random::RandomHemisphereVector( const vec3 &inDir )
{
	float f1 = fRandom(), f2 = fRandom();

	// find an orthonormal basis
	vec3 up = inDir;
	up.Normalize();
	vec3 u = ( up.X() > 0.99 || up.X() < -0.99 ? vec3::YAxis() : vec3::XAxis() );
	vec3 v = up.Cross( u );
	u = v.Cross( up );

	return (sqrt(f1)*cos(2.0f*M_PI*f2))*u + (sqrt(f1)*sin(2.0f*M_PI*f2))*v + (sqrt(1-f1))*up;

}


// used to record the random parameter for 
vec3 Random::RandomHemisphereVector( const vec3 &inDir, float &sample1, float &sample2 )
{
    float f1 = fRandom(), f2 = fRandom();
    //record the sample
    sample1 = f1;
    sample2 = f2;
    // find an orthonormal basis
    vec3 up = inDir;
    up.Normalize();
    vec3 u = ( up.X() > 0.99 || up.X() < -0.99 ? vec3::YAxis() : vec3::XAxis() );
    vec3 v = up.Cross( u );
    u = v.Cross( up );

    return (sqrt(f1)*cos(2.0f*M_PI*f2))*u + (sqrt(f1)*sin(2.0f*M_PI*f2))*v + (sqrt(1-f1))*up;

}

vec3 Random::RandomCosineWeightedHemisphereVector( const vec3 &inDir )
{
	float f1 = fRandom(), f2 = fRandom();
	float phi = M_PI*f1;
	float theta = 2.0f*M_PI*f2;

	// find an orthonormal basis
	vec3 up = inDir;
	up.Normalize();
	vec3 u = ( up.X() > 0.99 || up.X() < -0.99 ? vec3::YAxis() : vec3::XAxis() );
	vec3 v = up.Cross( u );
	u = v.Cross( up );

	return sin(phi)*cos(theta)*u + sin(phi)*sin(theta)*v + cos(phi)*up;
}



