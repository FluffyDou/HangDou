/******************************************************************/
/* igluHalton2D.h                                                 */
/* -----------------------                                        */
/*                                                                */
/* The file defines a generator for a Halton sequence of quasi    */
/*    random numbers on 2D surfaces.  Requires two prime bases    */
/*    that should also be mutually prime.                         */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUHALTON2D_H
#define IGLUHALTON2D_H

#include "igluHalton1D.h"
#include "igluSampler2D.h"

namespace iglu {

class IGLUHalton2D : public IGLUSampler2D
{
public:
	IGLUHalton2D( uint primeBase1=2, uint primeBase2=3 );
	~IGLUHalton2D() {}

	// The Halton sequence starts from sampIdx = 1
	//    (sampIdx=0 is NOT part of the sequence, though it will give a 
	//     valid sample [at 0,0] if you accidentally use it)
	virtual vec2 Sample( unsigned int sampIdx = 1 );

	// A pointer to a IGLUHalton2D could have type IGLUHalton2D::Ptr
	typedef IGLUHalton2D *Ptr;

private:
	IGLUHalton1D dim1, dim2;
};


// End iglu namespace
}


#endif

