/******************************************************************/
/* vec4.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 4D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_4D_H
#define IGLU_VECTOR_4D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class vec2;
class vec3;
class igluMatrix4x4;

class vec4 {
    float d[4];
    friend class vec2; 
	friend class vec3;
	friend class igluMatrix4x4;
public:
	// Vector constructors
	inline vec4(float allVals);                         // from a single float
	inline vec4(float data[4]);							// from a float array
    inline vec4(float x, float y, float z, float w);	// from individual floats
    inline vec4(const vec4& v);							// copy constructor from another vector
	inline vec4(const vec3& v, float);  				// copy constructor from another vector
	inline vec4(const vec2& v1, const vec2& v2);  		// copy constructor from another vector
	inline vec4(float, const vec3& v);  				// copy constructor from another vector
	inline vec4();										// default constructor

	// Operations for computing length and normalizing
	inline float Normalize();
	inline vec4 vNormalize() const;
    inline float Length( void ) const;
	inline float LengthSqr( void ) const;

	// Assignment constructors
    inline vec4& operator=(const vec4& v);

	// Boolean comparisons
    inline bool operator==(const vec4& v) const;
	inline bool operator != (const vec4& v) const;
    
	// Mathematical operators
	inline vec4 operator*(float s) const;
	friend inline vec4 operator*(float s, const vec4& v);
    inline vec4 operator*(const vec4& v) const;
    inline vec4 operator/(const vec4& v) const;
    inline vec4& operator*=(float s);
    inline vec4 operator+(const vec4& v) const;
    inline vec4& operator+=(const vec4& v);
    inline vec4 operator-() const;
    inline vec4 operator-(const vec4& v) const;

	// Dot and cross products
    inline vec4 Cross(const vec4& v) const;
    inline float Dot(const vec4& v) const;

	// Accessor functions
	inline float X() const						{ return d[0]; }
	inline float Y() const						{ return d[1]; }
	inline float Z() const						{ return d[2]; }
	inline float W() const						{ return d[3]; }
	inline float GetElement( int i ) const		{ return d[i]; }
	inline float *GetDataPtr()					{ return d; }
	inline const float *GetConstDataPtr() const	{ return d; }

	// Complex accessor functions
	inline vec3 xyz() const                     { return vec3( d[0], d[1], d[2] ); }
	inline vec3 xzy() const                     { return vec3( d[0], d[2], d[0] ); }
	inline vec3 zxy() const                     { return vec3( d[2], d[0], d[1] ); }
	inline vec2 xy() const                      { return vec2( d[0], d[1] ); }
	inline vec2 yz() const                      { return vec2( d[1], d[2] ); }
	inline vec2 zw() const                      { return vec2( d[2], d[3] ); }
	inline vec2 xz() const                      { return vec2( d[0], d[2] ); }
	inline vec2 xw() const                      { return vec2( d[0], d[3] ); }
	inline vec2 yw() const                      { return vec2( d[1], d[3] ); }

	// Static methods.  Useful for defining commonly used vectors
	static vec4 Zero( void )  { return vec4(0,0,0,0); }
	static vec4 One( void )   { return vec4(1,1,1,1); }
	static vec4 XAxis( void ) { return vec4(1,0,0,0); }
	static vec4 YAxis( void ) { return vec4(0,1,0,0); }
	static vec4 ZAxis( void ) { return vec4(0,0,1,0); }
	static vec4 WAxis( void ) { return vec4(0,0,0,1); }

	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend vec4 Min(const vec4& v1, const vec4& v2);
    friend vec4 Max(const vec4& v1, const vec4& v2);

	// A pretty-syntax normalization method
	friend vec3 Normalize( const vec3& v );

	// Computes a scalar triple product
	friend float ScalarTripleProduct(const vec4& v1, const vec4& v2, const vec4 &v3);

	// Returns the maximum (or minimum component)
	inline float MaxComponent( void ) const;
	inline float MinComponent( void ) const;

	inline void Print( void ) const { printf( "vec4: %f %f %f %f\n", d[0], d[1], d[2], d[3] ); }

	// A pointer to a vec4 could have type vec4::Ptr
    typedef vec4 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef vec4 fvec4;
typedef vec4 float4;

inline vec4::vec4(float allVals)
{
	d[0]=d[1]=d[2]=d[3]=allVals;
}

inline vec4::vec4( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline vec4::vec4(float data[4])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=data[2];
	d[3]=data[3];
}

inline vec4::vec4(float x, float y, float z, float w) 
{
    d[0]=x;
    d[1]=y;
    d[2]=z;
	d[3]=w;
}

inline vec4::vec4(const vec4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=v.d[3];
}

inline vec4::vec4(const vec3& v, float w)
{
	d[0]=v.d[0];
	d[1]=v.d[1];
	d[2]=v.d[2];
	d[3]=w;
}

inline vec4::vec4(const vec2& v1, const vec2& v2)
{
	d[0]=v1.d[0];
	d[1]=v1.d[1];
	d[2]=v2.d[0];
	d[3]=v2.d[1];
}

inline vec4::vec4(float x, const vec3& v)
{
	d[0]=x;
	d[1]=v.d[0];
	d[2]=v.d[1];
	d[3]=v.d[2];
}

inline float vec4::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
}

inline float vec4::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3];
}

inline vec4& vec4::operator=(const vec4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=v.d[3];
    return *this;
}

inline vec4 vec4::operator*(float s) const 
{
    return vec4(d[0]*s, d[1]*s, d[2]*s, d[3]*s);
}

inline vec4 operator*(float s, const vec4& v) 
{
	return vec4(v.d[0]*s, v.d[1]*s, v.d[2]*s, v.d[3]*s);
}


inline vec4 vec4::operator*(const vec4& v) const 
{
    return vec4(d[0]*v.d[0], d[1]*v.d[1], d[2]*v.d[2], d[3]*v.d[3]);
}

inline vec4 vec4::operator/(const vec4& v) const 
{
    return vec4(d[0]/v.d[0], d[1]/v.d[1], d[2]/v.d[2], d[3]/v.d[3]);
}

inline vec4 vec4::operator+(const vec4& v) const 
{
    return vec4(d[0]+v.d[0], d[1]+v.d[1], d[2]+v.d[2], d[3]+v.d[3]);
}

inline vec4& vec4::operator+=(const vec4& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    d[2]=d[2]+v.d[2];
	d[3]=d[3]+v.d[3];
    return *this;
}

inline vec4& vec4::operator*=(float s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    d[2]=d[2]*s;
	d[3]=d[3]*s;
    return *this;
}

inline vec4 vec4::operator-() const 
{
    return vec4(-d[0], -d[1], -d[2], -d[3]);
}

inline vec4 vec4::operator-(const vec4& v) const 
{
    return vec4(d[0]-v.d[0], d[1]-v.d[1], d[2]-v.d[2], d[3]-v.d[3]);
}

inline float vec4::Normalize() 
{
    float l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    d[2]=d[2]/l;
	d[3]=d[3]/l;
    return l;
}

inline vec4 vec4::vNormalize() const
{
	float l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
	return vec4( d[0]/l, d[1]/l, d[2]/l, d[3]/l );
}

inline vec4 vec4::Cross(const vec4& v) const 
{
    return vec4(d[1]*v.d[2]-d[2]*v.d[1],
    	        d[2]*v.d[0]-d[0]*v.d[2],
    	        d[0]*v.d[1]-d[1]*v.d[0],
				0.0);
}

inline float vec4::Dot(const vec4& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1]+d[2]*v.d[2]+d[3]*v.d[3];
}


inline bool vec4::operator != (const vec4& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2] || d[3] != v.d[3];
}

inline bool vec4::operator == (const vec4& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2] && d[3] == v.d[3];
}

inline float vec4::MinComponent( void ) const 
{
    float temp = (d[1] < d[0]? d[1] : d[0]); 
	temp = (d[2] < temp ? d[2] : temp );
	return (d[3] < temp ? d[3] : temp); 
}

inline float vec4::MaxComponent( void ) const 
{
	float temp = (d[1] > d[0]? d[1] : d[0]); 
	temp = (d[2] > temp ? d[2] : temp );
	return (d[3] > temp ? d[3] : temp); 
}


inline vec4 Min(const vec4& v1, const vec4& v2)
{
    return vec4( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]), min(v1.d[2], v2.d[2]), min(v1.d[3], v2.d[3]) );
}

inline vec4 Max(const vec4& v1, const vec4& v2)
{
    return vec4( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]), max(v1.d[2], v2.d[2]), max(v1.d[3], v2.d[3]) );
}

inline vec4 Normalize( const vec4& v )
{
	return v.vNormalize();
}

inline float ScalarTripleProduct(const vec4& v1, const vec4& v2, const vec4 &v3)
{
	return v1.Dot( v2.Cross( v3 ) );
}



// End iglu namespace
}

#endif

