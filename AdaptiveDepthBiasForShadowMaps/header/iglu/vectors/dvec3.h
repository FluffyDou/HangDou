/******************************************************************/
/* dvec3.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 3D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_DOUBLEVECTOR_3D_H
#define IGLU_DOUBLEVECTOR_3D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class dvec2;
class dvec4;

class dvec3 {
    double d[4]; 
	friend class dvec2;
	friend class dvec4;
public:
	// Vector constructors
	inline dvec3(double allVals);
	inline dvec3(double data[3]);                       // from a float array
    inline dvec3(double x, double y, double z);         // from individual floats
    inline dvec3(const dvec3& v);					    // copy constructor from another vector
	inline dvec3();                                     // default constructor

	// Operations for computing length and normalizing
	inline double Normalize();
	inline dvec3 vNormalize();
    inline double Length( void ) const;
	inline double LengthSqr( void ) const;

	// Assignment constructors
    inline dvec3& operator=(const dvec3& v);

	// Boolean comparisons
    inline bool operator==(const dvec3& v) const;
	inline bool operator != (const dvec3& v) const;
    
	// Mathematical operators
	inline dvec3 operator*(double s) const;
	friend inline dvec3 operator*(double s, const dvec3& v);
    inline dvec3 operator*(const dvec3& v) const;
    inline dvec3 operator/(const dvec3& v) const;
    inline dvec3& operator*=(double s);
    inline dvec3 operator+(const dvec3& v) const;
    inline dvec3& operator+=(const dvec3& v);
    inline dvec3 operator-() const;
    inline dvec3 operator-(const dvec3& v) const;

	// Dot and cross products
    inline dvec3 Cross(const dvec3& v) const;
    inline double Dot(const dvec3& v) const;

	// Accessor functions
	inline double X() const						 { return d[0]; }
	inline double Y() const						 { return d[1]; }
	inline double Z() const						 { return d[2]; }
	inline double GetElement( int i ) const		 { return d[i]; }
	inline double *GetDataPtr()					 { return d; }
	inline const double *GetConstDataPtr() const { return d; }

	// Static methods.  Useful for defining commonly used vectors
	static dvec3 Zero( void )  { return dvec3(0,0,0); }
	static dvec3 One( void )   { return dvec3(1,1,1); }
	static dvec3 XAxis( void ) { return dvec3(1,0,0); }
	static dvec3 YAxis( void ) { return dvec3(0,1,0); }
	static dvec3 ZAxis( void ) { return dvec3(0,0,1); }

	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend dvec3 Min(const dvec3& v1, const dvec3& v2);
    friend dvec3 Max(const dvec3& v1, const dvec3& v2);

	// Computes a scalar triple product
	friend double ScalarTripleProduct(const dvec3& v1, const dvec3& v2, const dvec3 &v3);

	// Returns the maximum (or minimum component)
	inline double MaxComponent( void ) const;
	inline double MinComponent( void ) const;

	inline void Print( void ) const { printf( "dvec3: %f %f %f\n", d[0], d[1], d[2] ); }

	// A pointer to a dvec3 could have type dvec3::Ptr
    typedef dvec3 *Ptr;
};

// Perhaps you want a consistent naming scheme among all vector types?
typedef dvec3 double3;


inline dvec3::dvec3(double allVals)
{
	d[0]=d[1]=d[2]=allVals;
	d[3]=0;
}

inline dvec3::dvec3( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline dvec3::dvec3(double data[3])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=data[2];
	d[3]=0;
}

inline dvec3::dvec3(double x, double y, double z) 
{
    d[0]=x;
    d[1]=y;
    d[2]=z;
	d[3]=0;
}

inline dvec3::dvec3(const dvec3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=0;
}

inline double dvec3::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
}

inline double dvec3::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
}

inline dvec3& dvec3::operator=(const dvec3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
    return *this;
}

inline dvec3 dvec3::operator*(double s) const 
{
    return dvec3(d[0]*s, d[1]*s, d[2]*s);
}

inline dvec3 operator*(double s, const dvec3& v) 
{
	return dvec3(v.d[0]*s, v.d[1]*s, v.d[2]*s);
}


inline dvec3 dvec3::operator*(const dvec3& v) const 
{
    return dvec3(d[0]*v.d[0], d[1]*v.d[1], d[2]*v.d[2]);
}

inline dvec3 dvec3::operator/(const dvec3& v) const 
{
    return dvec3(d[0]/v.d[0], d[1]/v.d[1], d[2]/v.d[2]);
}

inline dvec3 dvec3::operator+(const dvec3& v) const 
{
    return dvec3(d[0]+v.d[0], d[1]+v.d[1], d[2]+v.d[2]);
}

inline dvec3& dvec3::operator+=(const dvec3& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    d[2]=d[2]+v.d[2];
    return *this;
}

inline dvec3& dvec3::operator*=(double s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    d[2]=d[2]*s;
    return *this;
}

inline dvec3 dvec3::operator-() const 
{
    return dvec3(-d[0], -d[1], -d[2]);
}

inline dvec3 dvec3::operator-(const dvec3& v) const 
{
    return dvec3(d[0]-v.d[0], d[1]-v.d[1], d[2]-v.d[2]);
}

inline double dvec3::Normalize() 
{
    double l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    d[2]=d[2]/l;
    return l;
}

inline dvec3 dvec3::vNormalize()
{
	double l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	return dvec3( d[0]/l, d[1]/l, d[2]/l );
}

inline dvec3 dvec3::Cross(const dvec3& v) const 
{
    return dvec3(d[1]*v.d[2]-d[2]*v.d[1],
    	        d[2]*v.d[0]-d[0]*v.d[2],
    	        d[0]*v.d[1]-d[1]*v.d[0]);
}

inline double dvec3::Dot(const dvec3& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1]+d[2]*v.d[2];
}

inline bool dvec3::operator != (const dvec3& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2];
}

inline bool dvec3::operator == (const dvec3& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2];
}

inline double dvec3::MinComponent( void ) const 
{
    double temp = (d[1] < d[0]? d[1] : d[0]); 
	return (d[2] < temp? d[2] : temp); 
}

inline double dvec3::MaxComponent( void ) const 
{
	double temp = (d[1] > d[0]? d[1] : d[0]); 
	return (d[2] > temp? d[2] : temp); 
}


inline dvec3 Min(const dvec3& v1, const dvec3& v2)
{
    return dvec3( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]), min(v1.d[2], v2.d[2]) );
}

inline dvec3 Max(const dvec3& v1, const dvec3& v2)
{
    return dvec3( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]), max(v1.d[2], v2.d[2]) );
}

inline double ScalarTripleProduct(const dvec3& v1, const dvec3& v2, const dvec3 &v3)
{
	return v1.Dot( v2.Cross( v3 ) );
}


// End of namespace iglu
}



#endif

