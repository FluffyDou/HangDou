/******************************************************************/
/* dvec2.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_DOUBLEVECTOR_2D_H
#define IGLU_DOUBLEVECTOR_2D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class dvec3;
class dvec4;

class dvec2 {
    double d[4]; 
	friend class dvec3;
	friend class dvec4;
public:
	// Vector constructors
	inline dvec2(double allVals);
	inline dvec2(double data[2]);                      // from a float array
    inline dvec2(double x, double y);				   // from individual floats
    inline dvec2(const dvec2& v);		               // copy constructor from another vector
	inline dvec2();	                                   // default constructor

	// Operations for computing length and normalizing
	inline double Normalize();
	inline dvec2 vNormalize();
    inline double Length( void ) const;
	inline double LengthSqr( void ) const;

	// Assignment constructors
    inline dvec2& operator=(const dvec2& v);

	// Boolean comparisons
    inline bool operator==(const dvec2& v) const;
	inline bool operator != (const dvec2& v) const;
    
	// Mathematical operators
	inline dvec2 operator*(double s) const;
	friend inline dvec2 operator*(double s, const dvec2& v);
    inline dvec2 operator*(const dvec2& v) const;
    inline dvec2 operator/(const dvec2& v) const;
    inline dvec2& operator*=(double s);
    inline dvec2 operator+(const dvec2& v) const;
    inline dvec2& operator+=(const dvec2& v);
    inline dvec2 operator-() const;
    inline dvec2 operator-(const dvec2& v) const;

	// Dot and cross products
    inline double Dot(const dvec2& v) const;

	// Accessor functions
	inline double X() const						 { return d[0]; }
	inline double Y() const						 { return d[1]; }
	inline double GetElement( int i ) const		 { return d[i]; }
	inline double *GetDataPtr()					 { return d; }
	inline const double *GetConstDataPtr() const { return d; }

	// Static methods.  Useful for defining commonly used vectors
	static dvec2 Zero( void )				{ return dvec2(0,0); }
	static dvec2 One( void )				{ return dvec2(1,1); }
	static dvec2 XAxis( void )				{ return dvec2(1,0); }
	static dvec2 YAxis( void )				{ return dvec2(0,1); }


	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend dvec2 Min(const dvec2& v1, const dvec2& v2);
    friend dvec2 Max(const dvec2& v1, const dvec2& v2);

	// Returns the maximum (or minimum component)
	inline double MaxComponent( void ) const;
	inline double MinComponent( void ) const;

	inline void Print( void ) const { printf( "dvec2: %f %f\n", d[0], d[1] ); }

	// A pointer to a dvec2 could have type dvec2::Ptr
    typedef dvec2 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef dvec2 double2;


inline dvec2::dvec2(double allVals)
{
	d[0]=d[1]=allVals;
	d[2]=d[3]=0;
}

inline dvec2::dvec2( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline dvec2::dvec2(double data[2])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=0;
	d[3]=0;
}

inline dvec2::dvec2(double x, double y) 
{
    d[0]=x;
    d[1]=y;
    d[2]=0;
	d[3]=0;
}

inline dvec2::dvec2(const dvec2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=0;
	d[3]=0;
}

inline double dvec2::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]);
}

inline double dvec2::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1];
}

inline dvec2& dvec2::operator=(const dvec2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    return *this;
}

inline dvec2 dvec2::operator*(double s) const 
{
    return dvec2(d[0]*s, d[1]*s);
}

inline dvec2 operator*(double s, const dvec2& v) 
{
	return dvec2(v.d[0]*s, v.d[1]*s);
}


inline dvec2 dvec2::operator*(const dvec2& v) const 
{
    return dvec2(d[0]*v.d[0], d[1]*v.d[1]);
}

inline dvec2 dvec2::operator/(const dvec2& v) const 
{
    return dvec2(d[0]/v.d[0], d[1]/v.d[1]);
}

inline dvec2 dvec2::operator+(const dvec2& v) const 
{
    return dvec2(d[0]+v.d[0], d[1]+v.d[1]);
}

inline dvec2& dvec2::operator+=(const dvec2& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    return *this;
}

inline dvec2& dvec2::operator*=(double s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    return *this;
}

inline dvec2 dvec2::operator-() const 
{
    return dvec2(-d[0], -d[1]);
}

inline dvec2 dvec2::operator-(const dvec2& v) const 
{
    return dvec2(d[0]-v.d[0], d[1]-v.d[1]);
}

inline double dvec2::Normalize() 
{
    double l=sqrt(d[0]*d[0]+d[1]*d[1]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    return l;
}

inline dvec2 dvec2::vNormalize()
{
	double l=sqrt(d[0]*d[0]+d[1]*d[1]);
	return dvec2( d[0]/l, d[1]/l );
}

inline double dvec2::Dot(const dvec2& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1];
}

inline bool dvec2::operator != (const dvec2& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1];
}

inline bool dvec2::operator == (const dvec2& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1];
}

inline double dvec2::MinComponent( void ) const 
{
    return (d[1] < d[0]? d[1] : d[0]); 
}

inline double dvec2::MaxComponent( void ) const 
{
	return (d[1] > d[0]? d[1] : d[0]); 
}


inline dvec2 Min(const dvec2& v1, const dvec2& v2)
{
    return dvec2( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]) );
}

inline dvec2 Max(const dvec2& v1, const dvec2& v2)
{
    return dvec2( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]) );
}


// End of namespace iglu
}


#endif

