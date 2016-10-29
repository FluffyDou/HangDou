/******************************************************************/
/* uint2.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_INT_4D_H
#define IGLU_VECTOR_INT_4D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class int2;
class int3;

class int4 {
    int d[4]; 

public:
	// Vector constructors
	inline int4(int allVals);
	inline int4(int data[4]);                
    inline int4(int x, int y, int z, int w);	 
    inline int4(const int4& v);		               
	inline int4();	                                  

	// Assignment constructors
    inline int4& operator=(const int4& v);

	// Boolean comparisons
    inline bool operator==(const int4& v) const;
	inline bool operator != (const int4& v) const;
    
	// Accessor functions
	inline int X() const						{ return d[0]; }
	inline int Y() const						{ return d[1]; }
	inline int Z() const						{ return d[2]; }
	inline int W() const                        { return d[3]; }
	inline int GetElement( int i ) const		{ return d[i]; }
	inline int *GetDataPtr()					{ return d; }
	inline const int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "int4: %d %d %d %d\n", d[0], d[1], d[2], d[3] ); }

	// A pointer to a int4 could have type int4::Ptr
    typedef int4 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef int4 ivec4;

inline int4::int4(int allVals)
{
	d[0]=d[1]=d[2]=d[3]=allVals;
}

inline int4::int4( void )
{
    d[0]=d[1]=d[2]=d[3]=0;
}

inline int4::int4(int data[4])
{
    d[0]=data[0];
    d[1]=data[1];
	d[2]=data[2];
	d[3]=data[3];
}

inline int4::int4(int x, int y, int z, int w) 
{
    d[0]=x;
    d[1]=y;
	d[2]=z;
	d[3]=w;
}

inline int4::int4(const int4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
	d[3]=v.d[3];
}

inline int4& int4::operator=(const int4& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
	d[3]=v.d[3];
    return *this;
}

inline bool int4::operator != (const int4& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2] || d[3] != v.d[3];
}

inline bool int4::operator == (const int4& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2] && d[3] == v.d[3];
}


}

#endif

