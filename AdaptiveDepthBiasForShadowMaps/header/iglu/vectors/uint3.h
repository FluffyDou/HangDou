/******************************************************************/
/* uint2.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_UINT_3D_H
#define IGLU_VECTOR_UINT_3D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class uint2;
class uint4;

class uint3 {
    unsigned int d[3]; 

public:
	// Vector constructors
	inline uint3(unsigned int allVals);
	inline uint3(unsigned int data[3]);                
    inline uint3(unsigned int x, unsigned int y, unsigned int z);	 
    inline uint3(const uint3& v);		               
	inline uint3();	                                  

	// Assignment constructors
    inline uint3& operator=(const uint3& v);

	// Boolean comparisons
    inline bool operator==(const uint3& v) const;
	inline bool operator != (const uint3& v) const;
    
	// Accessor functions
	inline unsigned int X() const						{ return d[0]; }
	inline unsigned int Y() const						{ return d[1]; }
	inline unsigned int Z() const						{ return d[2]; }
	inline unsigned int GetElement( int i ) const		{ return d[i]; }
	inline unsigned int *GetDataPtr()					{ return d; }
	inline const unsigned int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "uint3: %d %d %d\n", d[0], d[1], d[2] ); }

	// A pointer to a uint3 could have type uint3::Ptr
    typedef uint3 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef uint3 uivec3;

inline uint3::uint3(unsigned int allVals)
{
	d[0]=d[1]=d[2]=allVals;
	d[3]=0;
}


inline uint3::uint3( void )
{
    d[0]=d[1]=d[2]=0;
}

inline uint3::uint3(unsigned int data[3])
{
    d[0]=data[0];
    d[1]=data[1];
	d[2]=data[2];
}

inline uint3::uint3(unsigned int x, unsigned int y, unsigned int z) 
{
    d[0]=x;
    d[1]=y;
	d[2]=z;
}

inline uint3::uint3(const uint3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
}

inline uint3& uint3::operator=(const uint3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
    return *this;
}

inline bool uint3::operator != (const uint3& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2];
}

inline bool uint3::operator == (const uint3& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2];
}


}

#endif

