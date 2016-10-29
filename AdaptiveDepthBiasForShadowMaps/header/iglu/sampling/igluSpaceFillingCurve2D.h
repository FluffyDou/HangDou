/******************************************************************/
/* igluSpaceFillingCurve2D.h                                      */
/* -----------------------                                        */
/*                                                                */
/* The file defines a base class for space-filling curves in 2D.  */
/*    (Basically these sorts of curves convert from 2D (x,y) to   */
/*    a single integer, and vice-versa)                           */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUSPACEFILLINGCURVE2D_H
#define IGLUSPACEFILLINGCURVE2D_H

#include "igluVector.h"

namespace iglu {

class IGLUSpaceFillingCurve2D
{
public:
	IGLUSpaceFillingCurve2D()  {}
	~IGLUSpaceFillingCurve2D() {}

	// Convert from (x,y) position to ordered position on the curve
	virtual int ToIndex( int2 xyPos ) = 0;
	
	// Convert from ordered position on the curve to (x,y) position
	virtual int2 ToXYCoord( int curveOrder ) = 0;

	// A pointer to a IGLUHalton2D could have type IGLUHalton2D::Ptr
	typedef IGLUSpaceFillingCurve2D *Ptr;
};


// End iglu namespace
}


#endif