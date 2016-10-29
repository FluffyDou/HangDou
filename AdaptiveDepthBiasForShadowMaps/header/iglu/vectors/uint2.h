/******************************************************************/
/* uint2.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_UINT_2D_H
#define IGLU_VECTOR_UINT_2D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class uint3;
class uint4;

class uint2 {
    unsigned int d[2]; 

public:
	// Vector constructors
	inline uint2(unsigned int allVals);
	inline uint2(unsigned int data[2]);                // from a float array
    inline uint2(unsigned int x, unsigned int y);	   // from individual floats
    inline uint2(const uint2& v);		               // copy constructor from another vector
	inline uint2();	                                   // default constructor

	// Assignment constructors
    inline uint2& operator=(const uint2& v);

	// Boolean comparisons
    inline bool operator==(const uint2& v) const;
	inline bool operator != (const uint2& v) const;
    
	// Accessor functions
	inline unsigned int X() const						{ return d[0]; }
	inline unsigned int Y() const						{ return d[1]; }
	inline unsigned int GetElement( int i ) const		{ return d[i]; }
	inline unsigned int *GetDataPtr()					{ return d; }
	inline const unsigned int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "uint2: %d %d\n", d[0], d[1] ); }

	// A pointer to a uint2 could have type uint2::Ptr
    typedef uint2 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef uint2 uivec2;

// On Linux/Unix systems uint is already typedef'd in "sys/types.h"
#ifndef __GNUC__
	typedef unsigned int uint;
#endif

typedef unsigned char  uchar;
typedef unsigned long  ulong;
typedef unsigned short ushort;


inline uint2::uint2(unsigned int allVals)
{
	d[0]=d[1]=allVals;
	d[2]=d[3]=0;
}

inline uint2::uint2( void )
{
    d[0]=d[1]=0;
}

inline uint2::uint2(unsigned int data[2])
{
    d[0]=data[0];
    d[1]=data[1];
}

inline uint2::uint2(unsigned int x, unsigned int y) 
{
    d[0]=x;
    d[1]=y;
}

inline uint2::uint2(const uint2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
}

inline uint2& uint2::operator=(const uint2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    return *this;
}

inline bool uint2::operator != (const uint2& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1];
}

inline bool uint2::operator == (const uint2& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1];
}


}


#endif

