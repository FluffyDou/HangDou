/******************************************************************/
/* igluMortonCurve2D.h                                            */
/* -----------------------                                        */
/*                                                                */
/* The file defines a class for converting from (x,y) to and from */
/*     a Morton ordering.  Requires the maximal ranges of x and y */
/*     to be powers of 2.                                         */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUMORTONCURVE2D_H
#define IGLUMORTONCURVE2D_H

#include "igluSpaceFillingCurve2D.h"

namespace iglu {

class IGLUMortonCurve2D : public IGLUSpaceFillingCurve2D
{
public:
	// Initialize a 2D Morton Curve indexer based on maximum x/y resolution (must be power of 2)
	IGLUMortonCurve2D( int xyResolution ) : IGLUSpaceFillingCurve2D(), m_resolution(xyResolution) {}
	~IGLUMortonCurve2D() {}

	// Convert from (x,y) position to ordered position on the curve
	virtual int ToIndex( int2 xyPos );
	
	// Convert from ordered position on the curve to (x,y) position
	virtual int2 ToXYCoord( int curveOrder );

	// A pointer to a IGLUHalton2D could have type IGLUHalton2D::Ptr
	typedef IGLUMortonCurve2D *Ptr;

private:
	int m_resolution;  // Max values of (x,y) in ToIndex() and ToXYCoord();

	// Tables for table-lookup dilation
	static unsigned short int dilateTable[256];
	static unsigned char      undilateTable[256];

	// See 'Converting to and from Dilated Integers' in IEEE Trans. Comput.
	//inline uint Dilate( ushort x )   { return dilateTable[ 0xFF & x ] | (dilateTable[ (0xFFFF & x) >> 8 ] << 16); }
	//inline ushort Undilate( uint x ) { return undilateTable[ 0xFF & ((x>>7)|x) ] | (undilateTable[ 0xFF & (((x>>7)|x)>>16) ] << 8); }

	inline uint Dilate( ushort x )   { uint r=x; r=(r|(r<<8))&0x00FF00FF; r=(r|(r<<4))&0x0F0F0F0F; r=(r|(r<<2))&0x33333333; r=(r|(r<<1))&0x55555555; return r; }
	inline ushort Undilate( uint x ) { x=(x*3)&0x66666666; x=(x*5)&0x78787878; x=(x*17)&0x7F807F80;x=(x*257)&0x7FFF8000; return ushort(x>>15); }
};


// End iglu namespace
}


#endif