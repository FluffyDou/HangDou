/******************************************************************/
/* int2.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 2D vector class that implments most   */
/*     operations you need to perform on vectors.                 */
/*                                                                */
/* Chris Wyman (06/08/2010)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_INT_2D_H
#define IGLU_VECTOR_INT_2D_H

#include <cmath>
#include <cstdio>
#include "igluMath.h"

namespace iglu {

class int3;
class int4;

class int2 {
    int d[2]; 

public:
	// Vector constructors
	inline int2(int allVals);
	inline int2(int data[2]);                
    inline int2(int x, int y);	 
    inline int2(const int2& v);		               
	inline int2();	                                  

	// Assignment constructors
    inline int2& operator=(const int2& v);

	// Boolean comparisons
    inline bool operator==(const int2& v) const;
	inline bool operator != (const int2& v) const;
    
	// Accessor functions
	inline int X() const						{ return d[0]; }
	inline int Y() const						{ return d[1]; }
	inline int GetElement( int i ) const		{ return d[i]; }
	inline int *GetDataPtr()					{ return d; }
	inline const int *GetConstDataPtr() const	{ return d; }

	inline void Print( void ) const { printf( "int2: %d %d\n", d[0], d[1] ); }

	// A pointer to a int2 could have type int2::Ptr
    typedef int2 *Ptr;
};


// Perhaps you want a consistent naming scheme among all vector types?
typedef int2 ivec2;

inline int2::int2(int allVals)
{
	d[0]=d[1]=allVals;
	d[2]=d[3]=0;
}

inline int2::int2( void )
{
    d[0]=d[1]=0;
}

inline int2::int2(int data[2])
{
    d[0]=data[0];
    d[1]=data[1];
}

inline int2::int2(int x, int y) 
{
    d[0]=x;
    d[1]=y;
}

inline int2::int2(const int2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
}

inline int2& int2::operator=(const int2& v) 
{
    d[0]=v.d[0];
    d[1]=v.d[1];
    return *this;
}

inline bool int2::operator != (const int2& v) const {
    return d[0] != v.d[0] || d[1] != v.d[1];
}

inline bool int2::operator == (const int2& v) const {
    return d[0] == v.d[0] && d[1] == v.d[1];
}


}


#endif

