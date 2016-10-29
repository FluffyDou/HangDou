/******************************************************************/
/* uint2.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_UINT_4D_H
#define IGLU_VECTOR_UINT_4D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class uint2;
class uint3;

class uint4 {
    unsigned int d[4]; 

public:
	// Vector constructors
	inline uint4(unsigned int allVals);
	inline uint4(unsigned int data[4]);                
    inline uint4(unsigned int x, unsigned int y, unsigned int z, unsigned int w);	 
    inline uint4(const uint4& v);		               
	inline uint4();	                                  

	// Assignment constructors
    inline uint4& operator=(const uint4& v);

	// Boolean comparisons
    inline bool operator==(const uint4& v) const;
	inline bool operator != (const uint4& v) const;
    
	// Accessor functions
	inline unsigned int X() const						{ return d[0]; }
	inline unsigned int Y() const						{ return d[1]; }
	inline unsigned int Z() const						{ return d[2]; }
	inline unsigned int W() const                       { return d[3]; }
	inline unsigned int GetElement( int i ) const		{ return d[i]; }
	inline unsigned int *GetDataPtr()					{ return d; }
	inline const unsigned int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "uint4: %d %d %d %d\n", d[0], d[1], d[2], d[3] ); }

	// A pointer to a uint4 could have type uint4::Ptr
    typedef uint4 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef uint4 uivec4;

inline uint4::uint4(unsigned int allVals)
{
	d[0]=d[1]=d[2]=d[3]=allVals;
}

inline uint4::uint4( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline uint4::uint4(unsigned int data[4])
{
    d[0]=data[0];
    d[1]=data[1];
	d[2]=data[2];
	d[3]=data[3];
}

inline uint4::uint4(unsigned int x, unsigned int y, unsigned int z, unsigned int w) 
{
    d[0]=x;
    d[1]=y;
	d[2]=z;
	d[3]=w;
}

inline uint4::uint4(const uint4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
	d[3]=v.d[3];
}

inline uint4& uint4::operator=(const uint4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
	d[3]=v.d[3];
    return *this;
}

inline bool uint4::operator != (const uint4& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2] || d[3] != v.d[3];
}

inline bool uint4::operator == (const uint4& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2] && d[3] == v.d[3];
}


}

#endif

