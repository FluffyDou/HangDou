/******************************************************************/
/* igluOrthoNormalBasis.h                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines an ortho normal basis class.                  */
/*                                                                */
/* This class takes in a vector and creates two mutually-         */
/*     perpendicular vectors that are both, also, perpendicular   */
/*     to the input vector.  In other words, it creates an ortho- */
/*     normal basis from the the input, and defines (in the new   */
/*     space) the input vector as the z-vector for the purposes   */
/*     of this class API.                                         */
/*                                                                */
/* Please note that this ortho normal basis is NOT unique, and if */  
/*     you pass in slightly varying inputs to the constructor,    */
/*     you may get wildly varying x and y axes.  However, often   */
/*     just having any basis is all that you need.                */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLU_ONB_H
#define IGLU_ONB_H

namespace iglu {

class vec3;

class IGLUOrthoNormalBasis
{
public:
	IGLUOrthoNormalBasis( const vec3 &zAxis );
	~IGLUOrthoNormalBasis() {}

	// This works with points or vectors in the same space as the
	//   original vector (zAxis) passed into the constructor.
	vec3 ToBasis( const vec3 &p );

	// This is not the way I envision using this.  But in case we need
	//     access to the individual basis vectors in the original space
	const vec3 &GetXAxis( void ) const          { return m_xAxis; }
	const vec3 &GetYAxis( void ) const          { return m_yAxis; }
	const vec3 &GetZAxis( void ) const          { return m_zAxis; }

private:
	// The vectors for the orthonormal basis.
	vec3 m_xAxis, m_yAxis, m_zAxis;
};


// End iglu namespace
}


#endif

