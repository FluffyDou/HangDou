/******************************************************************/
/* vec3.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 3D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_3D_H
#define IGLU_VECTOR_3D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class vec2;
class vec4;

class vec3 {
    float d[4]; 
	friend class vec2;
	friend class vec4;
public:
	// Vector constructors
	inline vec3(float allVals);                        // from a single float
	inline vec3(float data[3]);                        // from a float array
    inline vec3(float x, float y, float z);            // from individual floats
    inline vec3(const vec3& v);						    // copy constructor from another vector
	inline vec3();                                     // default constructor

	// Operations for computing length and normalizing
	inline float Normalize();
	inline vec3 vNormalize() const;
    inline float Length( void ) const;
	inline float LengthSqr( void ) const;

	// Assignment constructors
    inline vec3& operator=(const vec3& v);

	// Boolean comparisons
    inline bool operator==(const vec3& v) const;
	inline bool operator != (const vec3& v) const;
    
	// Mathematical operators
	inline vec3 operator*(float s) const;
	friend inline vec3 operator*(float s, const vec3& v);
    inline vec3 operator*(const vec3& v) const;
    inline vec3 operator/(const vec3& v) const;
    inline vec3& operator*=(float s);
    inline vec3 operator+(const vec3& v) const;
    inline vec3& operator+=(const vec3& v);
    inline vec3 operator-() const;
    inline vec3 operator-(const vec3& v) const;

	// Dot and cross products
    inline vec3 Cross(const vec3& v) const;
    inline float Dot(const vec3& v) const;

	// Accessor functions
	inline float X() const						{ return d[0]; }
	inline float Y() const						{ return d[1]; }
	inline float Z() const						{ return d[2]; }
	inline float GetElement( int i ) const		{ return d[i]; }
	inline float *GetDataPtr()					{ return d; }
	inline const float *GetConstDataPtr() const	{ return d; }

	// Static methods.  Useful for defining commonly used vectors
	static vec3 Zero( void )  { return vec3(0,0,0); }
	static vec3 One( void )   { return vec3(1,1,1); }
	static vec3 XAxis( void ) { return vec3(1,0,0); }
	static vec3 YAxis( void ) { return vec3(0,1,0); }
	static vec3 ZAxis( void ) { return vec3(0,0,1); }

	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend vec3 Min(const vec3& v1, const vec3& v2);
    friend vec3 Max(const vec3& v1, const vec3& v2);

	// A pretty-syntax normalization method
	friend vec3 Normalize( const vec3& v );

	// Computes a scalar triple product
	friend float ScalarTripleProduct(const vec3& v1, const vec3& v2, const vec3 &v3);

	// Returns the maximum (or minimum component)
	inline float MaxComponent( void ) const;
	inline float MinComponent( void ) const;

	inline void Print( void ) const { printf( "vec3: %f %f %f\n", d[0], d[1], d[2] ); }

	// A pointer to a vec3 could have type vec3::Ptr
    typedef vec3 *Ptr;
};

// Perhaps you want a consistent naming scheme among all vector types?
typedef vec3 fvec3;
typedef vec3 float3;

inline vec3::vec3(float allVals)
{
	d[0]=d[1]=d[2]=allVals;
	d[3]=0;
}

inline vec3::vec3( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline vec3::vec3(float data[3])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=data[2];
	d[3]=0;
}

inline vec3::vec3(float x, float y, float z) 
{
    d[0]=x;
    d[1]=y;
    d[2]=z;
	d[3]=0;
}

inline vec3::vec3(const vec3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=0;
}

inline float vec3::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
}

inline float vec3::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
}

inline vec3& vec3::operator=(const vec3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
    return *this;
}

inline vec3 vec3::operator*(float s) const 
{
    return vec3(d[0]*s, d[1]*s, d[2]*s);
}

inline vec3 operator*(float s, const vec3& v) 
{
	return vec3(v.d[0]*s, v.d[1]*s, v.d[2]*s);
}


inline vec3 vec3::operator*(const vec3& v) const 
{
    return vec3(d[0]*v.d[0], d[1]*v.d[1], d[2]*v.d[2]);
}

inline vec3 vec3::operator/(const vec3& v) const 
{
    return vec3(d[0]/v.d[0], d[1]/v.d[1], d[2]/v.d[2]);
}

inline vec3 vec3::operator+(const vec3& v) const 
{
    return vec3(d[0]+v.d[0], d[1]+v.d[1], d[2]+v.d[2]);
}

inline vec3& vec3::operator+=(const vec3& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    d[2]=d[2]+v.d[2];
    return *this;
}

inline vec3& vec3::operator*=(float s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    d[2]=d[2]*s;
    return *this;
}

inline vec3 vec3::operator-() const 
{
    return vec3(-d[0], -d[1], -d[2]);
}

inline vec3 vec3::operator-(const vec3& v) const 
{
    return vec3(d[0]-v.d[0], d[1]-v.d[1], d[2]-v.d[2]);
}

inline float vec3::Normalize() 
{
    float l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    d[2]=d[2]/l;
    return l;
}

inline vec3 vec3::vNormalize() const
{
	float l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	return vec3( d[0]/l, d[1]/l, d[2]/l );
}

inline vec3 vec3::Cross(const vec3& v) const 
{
    return vec3(d[1]*v.d[2]-d[2]*v.d[1],
    	        d[2]*v.d[0]-d[0]*v.d[2],
    	        d[0]*v.d[1]-d[1]*v.d[0]);
}

inline float vec3::Dot(const vec3& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1]+d[2]*v.d[2];
}

inline bool vec3::operator != (const vec3& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2];
}

inline bool vec3::operator == (const vec3& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2];
}

inline float vec3::MinComponent( void ) const 
{
    float temp = (d[1] < d[0]? d[1] : d[0]); 
	return (d[2] < temp? d[2] : temp); 
}

inline float vec3::MaxComponent( void ) const 
{
	float temp = (d[1] > d[0]? d[1] : d[0]); 
	return (d[2] > temp? d[2] : temp); 
}


inline vec3 Min(const vec3& v1, const vec3& v2)
{
    return vec3( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]), min(v1.d[2], v2.d[2]) );
}

inline vec3 Max(const vec3& v1, const vec3& v2)
{
    return vec3( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]), max(v1.d[2], v2.d[2]) );
}

inline vec3 Normalize( const vec3& v )
{
	return v.vNormalize();
}

inline float ScalarTripleProduct(const vec3& v1, const vec3& v2, const vec3 &v3)
{
	return v1.Dot( v2.Cross( v3 ) );
}



// End iglu namespace
}

#endif

