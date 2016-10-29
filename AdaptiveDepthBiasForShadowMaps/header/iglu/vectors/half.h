/******************************************************************/
/* half.h                                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic 16-bit half float class               */
/*     (Currently extremely limited...  And not efficient.)       */
/*                                                                */
/* Chris Wyman (04/03/2013)                                       */
/******************************************************************/

#ifndef IGLU_HALF_H
#define IGLU_HALF_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include "igluMath.h"

namespace iglu {


class half {
    unsigned short data;
public:
	// Vector constructors
	       half( float v );                        // from a single float
	inline half( const half& v );                  // copy constructor from another half
	       half( void );                           // default constructor
	
	// Operations for pulling out other data representations
	       float          AsFloat( void ) const; 
	inline unsigned short AsShort( void ) const     { return data; }
	
	// Assignment constructors
    inline half& operator=(const half& v);

	// Boolean comparisons
    inline bool operator==(const half& v) const;
	inline bool operator != (const half& v) const;

	inline void Print( void ) const { printf( "half: %f\n", AsFloat() ); }

	// A pointer to a half could have type half::Ptr
    typedef half *Ptr;
};


inline half::half(const half& v) 
{
    data = v.data;
}

inline half& half::operator=(const half& v) 
{
    data = v.data;
    return *this;
}

inline bool half::operator != (const half& v) const {
    return data != v.data;
}

inline bool half::operator == (const half& v) const {
    return data == v.data;
}



// End iglu namespace
}

#endif

