/******************************************************************/
/* igluMath.h                                                     */
/* -----------------------                                        */
/*                                                                */
/* Some basic, common math functions not defined elsewhere.       */
/*                                                                */
/* Chris Wyman (12/06/2009)                                       */
/******************************************************************/

#ifndef IGLU_MATH
#define IGLU_MATH

// All iglu functions and classes are in the "iglu" namespace
namespace iglu
{

// Simple minimum and maximum functions
inline int min( int x, int y )		        { return x < y ? x : y; }
inline int max( int x, int y )		        { return x > y ? x : y; }
inline float min( float x, float y )		{ return x < y ? x : y; }
inline float max( float x, float y )		{ return x > y ? x : y; }
inline double min( double x, double y )		{ return x < y ? x : y; }
inline double max( double x, double y )     { return x > y ? x : y; }

// Function for computping log-base 2
float log2( float x );

// Function for squaring a value
inline float sqr( float x )					{ return x*x; }


// Constants that are often left undefined in standard headers
const double IGLU_PI = 3.141592653589793238462643;
const double IGLU_SQRT2 = 1.4142135623730950488016887;


// End iglu namespace
}

#endif