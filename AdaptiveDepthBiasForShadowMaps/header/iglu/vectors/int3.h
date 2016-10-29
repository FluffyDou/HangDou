/******************************************************************/
/* int3.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_INT_3D_H
#define IGLU_VECTOR_INT_3D_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class int2;
class int3;

class int3 {
    int d[3]; 

public:
	// Vector constructors
	inline int3(int allVals);
	inline int3(int data[3]);                
    inline int3(int x, int y, int z);	 
    inline int3(const int3& v);		               
	inline int3();	                                  

	// Assignment constructors
    inline int3& operator=(const int3& v);

	// Boolean comparisons
    inline bool operator==(const int3& v) const;
	inline bool operator != (const int3& v) const;
    
	// Accessor functions
	inline int X() const						{ return d[0]; }
	inline int Y() const						{ return d[1]; }
	inline int Z() const						{ return d[2]; }
	inline int GetElement( int i ) const		{ return d[i]; }
	inline int *GetDataPtr()					{ return d; }
	inline const int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "int3: %d %d %d\n", d[0], d[1], d[2] ); }

	// A pointer to a int3 could have type int3::Ptr
    typedef int3 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef int3 ivec3;

inline int3::int3(int allVals)
{
	d[0]=d[1]=d[2]=allVals;
	d[3]=0;
}


inline int3::int3( void )
{
    d[0]=d[1]=d[2]=0;
}

inline int3::int3(int data[3])
{
    d[0]=data[0];
    d[1]=data[1];
	d[2]=data[2];
}

inline int3::int3(int x, int y, int z) 
{
    d[0]=x;
    d[1]=y;
	d[2]=z;
}

inline int3::int3(const int3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
}

inline int3& int3::operator=(const int3& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
	d[2]=v.d[2];
    return *this;
}

inline bool int3::operator != (const int3& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1] || d[2] != v.d[2];
}

inline bool int3::operator == (const int3& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1] && d[2] == v.d[2];
}



}

#endif

