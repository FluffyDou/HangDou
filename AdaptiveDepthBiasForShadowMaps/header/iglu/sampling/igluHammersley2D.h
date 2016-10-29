/******************************************************************/
/* igluHammersley2D.h                                             */
/* -----------------------                                        */
/*                                                                */
/* The file defines a generator for a Hammersley sequence of      */
/*    quasi random numbers on 2D surfaces.  The Hammersley        */
/*    requires the user to know how many sample will be used      */
/*    before sampling is started.  In the case of this class,     */
/*    this means the number of samples needs to be passed in as   */
/*    input.                                                      */
/*                                                                */
/* Do NOT try to "cheat" by inputting a number of samples bigger  */
/*    that the number you expect to use.  The quality of the      */
/*    samples depends on using exactly the number specified (so   */
/*    if you say "I'll use 32 samples" and then use 16, you will  */
/*    only sample half the unit square)!                          */
/*                                                                */
/* If you want to use an arbitrary number of samples, use the     */
/*    Halton sequence rather than the Hammersley!                 */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUHAMMERSLEY2D_H
#define IGLUHAMMERSLEY2D_H

#include "igluHalton1D.h"
#include "igluSampler2D.h"

namespace iglu {

class IGLUHammersley2D : public IGLUSampler2D
{
public:
	// Input the total number of samples you will use during sampling
	IGLUHammersley2D( uint numSamples );
	~IGLUHammersley2D() {}

	// The Hammersley sequence starts from sampIdx = 1
	//    (sampIdx=0 is NOT part of the sequence, though it will give a 
	//     valid sample [at 0,0] if you accidentally use it)
	virtual vec2 Sample( unsigned int sampIdx = 1 );

	// A pointer to a IGLUHammersley2D could have type IGLUHammersley2D::Ptr
	typedef IGLUHammersley2D *Ptr;

private:
	IGLUHalton1D m_dim2;
	uint m_numSamples;
};


// End iglu namespace
}


#endif

