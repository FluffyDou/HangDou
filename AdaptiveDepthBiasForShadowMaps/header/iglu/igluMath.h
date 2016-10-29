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

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

// All iglu functions and classes are in the "iglu" namespace
namespace iglu
{

// Simple minimum and maximum functions
inline int min( int x, int y )		                      { return x < y ? x : y; }
inline int max( int x, int y )		                      { return x > y ? x : y; }
inline int Min( int x, int y )		                      { return x < y ? x : y; }
inline int Max( int x, int y )		                      { return x > y ? x : y; }
inline int clamp( int x, int _min, int _max )             { return (x < _min) ? _min : ((x>_max) ? _max : x ); }
inline float min( float x, float y )		              { return x < y ? x : y; }
inline float max( float x, float y )		              { return x > y ? x : y; }
inline float Min( float x, float y )		              { return x < y ? x : y; }
inline float Max( float x, float y )		              { return x > y ? x : y; }
inline float clamp( float x, float _min, float _max )     { return (x < _min) ? _min : ((x>_max) ? _max : x ); }
inline double min( double x, double y )		              { return x < y ? x : y; }
inline double max( double x, double y )                   { return x > y ? x : y; }
inline double Min( double x, double y )		              { return x < y ? x : y; }
inline double Max( double x, double y )                   { return x > y ? x : y; }
inline double clamp( double x, double _min, double _max ) { return (x < _min) ? _min : ((x>_max) ? _max : x ); }

// Function for computping log-base 2
float log2( float x );

// Function for squaring a value
inline int sqr( int x )                                 { return x*x; }
inline float sqr( float x )                             { return x*x; }
inline double sqr( double x )                           { return x*x; }


// Constants that are often left undefined in standard headers
const double IGLU_PI       = 3.14159265358979323846264338327950288;
const double IGLU_PI_2     = 1.57079632679489661923132169163975144;
const double IGLU_PI_4     = 0.785398163397448309615660845819875721;
const double IGLU_1_PI     = 0.318309886183790671537767526745028724;
const double IGLU_2_PI     = 0.636619772367581343075535053490057448;
const double IGLU_2_SQRTPI = 1.12837916709551257389615890312154517;
const double IGLU_SQRT2    = 1.41421356237309504880168872420969808;
const double IGLU_SQRT1_2  = 0.707106781186547524400844362104849039;
const double IGLU_E        = 2.71828182845904523536028747135266250;
const double IGLU_LOG2E    = 1.44269504088896340735992468100189214;
const double IGLU_LOG10E   = 0.434294481903251827651128918916605082;
const double IGLU_LN2      = 0.693147180559945309417232121458176568;
const double IGLU_LN10     = 2.30258509299404568401799145468436421;



// End iglu namespace
}

#endif