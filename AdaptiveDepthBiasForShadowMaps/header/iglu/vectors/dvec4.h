/******************************************************************/
/* dvec4.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 4D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_DOUBLEVECTOR_4D_H
#define IGLU_DOUBLEVECTOR_4D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class dvec2;
class dvec3;

class dvec4 {
    double d[4]; 
    friend class dvec2; 
	friend class dvec3;
public:
	// Vector constructors
	inline dvec4(double allVals);
	inline dvec4(double data[4]);							// from a double array
    inline dvec4(double x, double y, double z, double w);	// from individual floats
    inline dvec4(const dvec4& v);							// copy constructor from another vector
	inline dvec4();									        // default constructor

	// Operations for computing length and normalizing
	inline double Normalize();
	inline dvec4 vNormalize();
    inline double Length( void ) const;
	inline double LengthSqr( void ) const;

	// Assignment constructors
    inline dvec4& operator=(const dvec4& v);

	// Boolean comparisons
    inline bool operator==(const dvec4& v) const;
	inline bool operator != (const dvec4& v) const;
    
	// Mathematical operators
	inline dvec4 operator*(double s) const;
	friend inline dvec4 operator*(double s, const dvec4& v);
    inline dvec4 operator*(const dvec4& v) const;
    inline dvec4 operator/(const dvec4& v) const;
    inline dvec4& operator*=(double s);
    inline dvec4 operator+(const dvec4& v) const;
    inline dvec4& operator+=(const dvec4& v);
    inline dvec4 operator-() const;
    inline dvec4 operator-(const dvec4& v) const;

	// Dot and cross products
    inline dvec4 Cross(const dvec4& v) const;
    inline double Dot(const dvec4& v) const;

	// Accessor functions
	inline double X() const						 { return d[0]; }
	inline double Y() const						 { return d[1]; }
	inline double Z() const						 { return d[2]; }
	inline double W() const						 { return d[3]; }
	inline double GetElement( int i ) const		 { return d[i]; }
	inline double *GetDataPtr()					 { return d; }
	inline const double *GetConstDataPtr() const { return d; }

	// Static methods.  Useful for defining commonly used vectors
	static dvec4 Zero( void )  { return dvec4(0,0,0,0); }
	static dvec4 One( void )   { return dvec4(1,1,1,1); }
	static dvec4 XAxis( void ) { return dvec4(1,0,0,0); }
	static dvec4 YAxis( void ) { return dvec4(0,1,0,0); }
	static dvec4 ZAxis( void ) { return dvec4(0,0,1,0); }
	static dvec4 WAxis( void ) { return dvec4(0,0,0,1); }

	// Performs component-wise max and min operations on vectors
	//    (NOTE: These are not particularly fast!)
	friend dvec4 Min(const dvec4& v1, const dvec4& v2);
    friend dvec4 Max(const dvec4& v1, const dvec4& v2);

	// Computes a scalar triple product
	friend double ScalarTripleProduct(const dvec4& v1, const dvec4& v2, const dvec4 &v3);

	// Returns the maximum (or minimum component)
	inline double MaxComponent( void ) const;
	inline double MinComponent( void ) const;

	inline void Print( void ) const { printf( "dvec4: %f %f %f %f\n", d[0], d[1], d[2], d[3] ); }

	// A pointer to a dvec4 could have type dvec4::Ptr
    typedef dvec4 *Ptr;
};

// Perhaps you want a consistent naming scheme among all vector types?
typedef dvec4 double4;


inline dvec4::dvec4(double allVals)
{
	d[0]=d[1]=d[2]=d[3]=allVals;
}

inline dvec4::dvec4( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline dvec4::dvec4(double data[4])
{
    d[0]=data[0];
    d[1]=data[1];
    d[2]=data[2];
	d[3]=data[3];
}

inline dvec4::dvec4(double x, double y, double z, double w) 
{
    d[0]=x;
    d[1]=y;
    d[2]=z;
	d[3]=w;
}

inline dvec4::dvec4(const dvec4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=v.d[3];
}

inline double dvec4::Length( void ) const 
{
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
}

inline double dvec4::LengthSqr( void ) const 
{
    return d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3];
}

inline dvec4& dvec4::operator=(const dvec4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    d[2]=v.d[2];
	d[3]=v.d[3];
    return *this;
}

inline dvec4 dvec4::operator*(double s) const 
{
    return dvec4(d[0]*s, d[1]*s, d[2]*s, d[3]*s);
}

inline dvec4 operator*(double s, const dvec4& v) 
{
	return dvec4(v.d[0]*s, v.d[1]*s, v.d[2]*s, v.d[3]*s);
}


inline dvec4 dvec4::operator*(const dvec4& v) const 
{
    return dvec4(d[0]*v.d[0], d[1]*v.d[1], d[2]*v.d[2], d[3]*v.d[3]);
}

inline dvec4 dvec4::operator/(const dvec4& v) const 
{
    return dvec4(d[0]/v.d[0], d[1]/v.d[1], d[2]/v.d[2], d[3]/v.d[3]);
}

inline dvec4 dvec4::operator+(const dvec4& v) const 
{
    return dvec4(d[0]+v.d[0], d[1]+v.d[1], d[2]+v.d[2], d[3]+v.d[3]);
}

inline dvec4& dvec4::operator+=(const dvec4& v) 
{
    d[0]=d[0]+v.d[0];
    d[1]=d[1]+v.d[1];
    d[2]=d[2]+v.d[2];
	d[3]=d[3]+v.d[3];
    return *this;
}

inline dvec4& dvec4::operator*=(double s) 
{
    d[0]=d[0]*s;
    d[1]=d[1]*s;
    d[2]=d[2]*s;
	d[3]=d[3]*s;
    return *this;
}

inline dvec4 dvec4::operator-() const 
{
    return dvec4(-d[0], -d[1], -d[2], -d[3]);
}

inline dvec4 dvec4::operator-(const dvec4& v) const 
{
    return dvec4(d[0]-v.d[0], d[1]-v.d[1], d[2]-v.d[2], d[3]-v.d[3]);
}

inline double dvec4::Normalize() 
{
    double l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
    d[0]=d[0]/l;
    d[1]=d[1]/l;
    d[2]=d[2]/l;
	d[3]=d[3]/l;
    return l;
}

inline dvec4 dvec4::vNormalize()
{
	double l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
	return dvec4( d[0]/l, d[1]/l, d[2]/l, d[3]/l );
}

inline dvec4 dvec4::Cross(const dvec4& v) const 
{
    return dvec4(d[1]*v.d[2]-d[2]*v.d[1],
    	        d[2]*v.d[0]-d[0]*v.d[2],
    	        d[0]*v.d[1]-d[1]*v.d[0],
				0.0);
}

inline double dvec4::Dot(const dvec4& v) const 
{
    return d[0]*v.d[0]+d[1]*v.d[1]+d[2]*v.d[2]+d[3]*v.d[3];
}


inline bool dvec4::operator != (const dvec4& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2] || d[3] != v.d[3];
}

inline bool dvec4::operator == (const dvec4& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2] && d[3] == v.d[3];
}

inline double dvec4::MinComponent( void ) const 
{
    double temp = (d[1] < d[0]? d[1] : d[0]); 
	temp = (d[2] < temp ? d[2] : temp );
	return (d[3] < temp ? d[3] : temp); 
}

inline double dvec4::MaxComponent( void ) const 
{
	double temp = (d[1] > d[0]? d[1] : d[0]); 
	temp = (d[2] > temp ? d[2] : temp );
	return (d[3] > temp ? d[3] : temp); 
}


inline dvec4 Min(const dvec4& v1, const dvec4& v2)
{
    return dvec4( min(v1.d[0], v2.d[0]), min(v1.d[1], v2.d[1]), min(v1.d[2], v2.d[2]), min(v1.d[3], v2.d[3]) );
}

inline dvec4 Max(const dvec4& v1, const dvec4& v2)
{
    return dvec4( max(v1.d[0], v2.d[0]), max(v1.d[1], v2.d[1]), max(v1.d[2], v2.d[2]), max(v1.d[3], v2.d[3]) );
}

inline double ScalarTripleProduct(const dvec4& v1, const dvec4& v2, const dvec4 &v3)
{
	return v1.Dot( v2.Cross( v3 ) );
}


// End of namespace iglu
}

#endif

