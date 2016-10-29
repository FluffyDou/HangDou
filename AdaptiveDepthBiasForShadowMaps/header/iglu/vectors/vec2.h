/******************************************************************/
/* vec2.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_2D_H
#define IGLU_VECTOR_2D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class vec3;
class vec4;

class vec2 {
    float d[4]; 
	friend class vec3;
	friend class vec4;
public:
	// Vector constructors
	inline vec2(float allVals);                        // from a single float
	inline vec2(float data[2]);                        // from a float array
    inline vec2(float x, float y);					   // from individual floats
    inline vec2(const vec2& v);		                   // copy constructor from another vector
	inline vec2();	                                   // default constructor
	
	// Operations for computing length and normalizing
	inline float Normalize();
	inline vec2 vNormalize() const;
    inline float Length( void ) const;
	inline float LengthSqr( void ) const;

	// Assignment constructors
    inline vec2& operator=(const vec2& v);

	// Boolean comparisons
    inline bool operator==(const vec2& v) const;
	inline bool operator != (const vec2& v) const;
    
	// Mathematical operators
	inline vec2 operator*(float s) const;
	friend inline vec2 operator*(float s, const vec2& v);
    inline vec2 operator*(const vec2& v) const;
    inline vec2 operator/(const vec2& v) const;
    inline vec2& operator*=(float s);
    inline vec2 operator+(const vec2& v) const;
    inline vec2& operator+=(const vec2& v);
    inline vec2 operator-() const;
    inline vec2 operator-(const vec2& v) const;

	// Dot and cross products
    inline float Dot(const vec2& v) const;

	// Accessor functions
	inline float X() const						{ return d[0]; }
	inline float Y() const						{ return d[1]; }
	inline float GetElement( int i ) const		{ return d[i]; }
	inline float *GetDataPtr()					{ return d; }
	inline const float *GetConstDataPtr() const	{ return d; }

	// Static methods.  Useful for defining commonly used vectors
	static vec2 Zero( void )				{ return vec2(0,0); }
	static vec2 One( void )					{ return vec2(1,1); }
	static vec2 XAxis( void )				{ return vec2(1,0); }
	static vec2 YAxis( void )				{ return vec2(0,1); }


	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend vec2 Min(const vec2& v1, const vec2& v2);
    friend vec2 Max(const vec2& v1, const vec2& v2);

	// Returns the maximum (or minimum component)
	inline float MaxComponent( void ) const;
	inline float MinComponent( void ) const;

	// A pretty-syntax normalization method
	friend vec3 Normalize( const vec3& v );

	inline void Print( void ) const { printf( "vec2: %f %f\n", d[0], d[1] ); }

	// A pointer to a vec2 could have type vec2::Ptr
    typedef vec2 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef vec2 fvec2;
typedef vec2 float2;

inline vec2::vec2(float allVals)
{
	d[0]=d[1]=allVals;
	d[2]=d[3]=0;
}

inline vec2::vec2( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline vec2::vec2(float data[2])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=0;
	d[3]=0;
}

inline vec2::vec2(float x, float y) 
{
    d[0]=x;
    d[1]=y;
    d[2]=0;
	d[3]=0;
}

inline vec2::vec2(const vec2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=0;
	d[3]=0;
}

inline float vec2::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]);
}

inline float vec2::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1];
}

inline vec2& vec2::operator=(const vec2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    return *this;
}

inline vec2 vec2::operator*(float s) const 
{
    return vec2(d[0]*s, d[1]*s);
}

inline vec2 operator*(float s, const vec2& v) 
{
	return vec2(v.d[0]*s, v.d[1]*s);
}


inline vec2 vec2::operator*(const vec2& v) const 
{
    return vec2(d[0]*v.d[0], d[1]*v.d[1]);
}

inline vec2 vec2::operator/(const vec2& v) const 
{
    return vec2(d[0]/v.d[0], d[1]/v.d[1]);
}

inline vec2 vec2::operator+(const vec2& v) const 
{
    return vec2(d[0]+v.d[0], d[1]+v.d[1]);
}

inline vec2& vec2::operator+=(const vec2& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    return *this;
}

inline vec2& vec2::operator*=(float s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    return *this;
}

inline vec2 vec2::operator-() const 
{
    return vec2(-d[0], -d[1]);
}

inline vec2 vec2::operator-(const vec2& v) const 
{
    return vec2(d[0]-v.d[0], d[1]-v.d[1]);
}

inline float vec2::Normalize() 
{
    float l=sqrt(d[0]*d[0]+d[1]*d[1]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    return l;
}

inline vec2 vec2::vNormalize() const
{
	float l=sqrt(d[0]*d[0]+d[1]*d[1]);
	return vec2( d[0]/l, d[1]/l );
}

inline float vec2::Dot(const vec2& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1];
}

inline bool vec2::operator != (const vec2& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1];
}

inline bool vec2::operator == (const vec2& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1];
}

inline float vec2::MinComponent( void ) const 
{
    return (d[1] < d[0]? d[1] : d[0]); 
}

inline float vec2::MaxComponent( void ) const 
{
	return (d[1] > d[0]? d[1] : d[0]); 
}


inline vec2 Min(const vec2& v1, const vec2& v2)
{
    return vec2( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]) );
}

inline vec2 Max(const vec2& v1, const vec2& v2)
{
    return vec2( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]) );
}

inline vec2 Normalize( const vec2& v )
{
	return v.vNormalize();
}

// End iglu namespace
}

#endif

