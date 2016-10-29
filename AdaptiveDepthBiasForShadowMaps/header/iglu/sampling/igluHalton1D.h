/******************************************************************/
/* igluHalton1D.h                                                 */
/* -----------------------                                        */
/*                                                                */
/* The file defines a generator for a Halton sequence of quasi    */
/*    random numbers, given a particular prime number as base.    */
/*    For simple uses, the base will be small prime numbers (e.g. */
/*    2 for one dimension, 3 for another dimension, 5 for the 3rd */
/*    dimension) that vary from one sampling dimension to another */
/*    To ensure good sampling in high dimensional spaces, make    */
/*    there are some de-corellation techniques needed once you    */
/*    get beyond a base of 7 or 11.  This class does not handle   */
/*    these.                                                      */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUHALTON1D_H
#define IGLUHALTON1D_H

namespace iglu {

class IGLUHalton1D
{
public:
	IGLUHalton1D( unsigned int primeBase=2 );
	~IGLUHalton1D() {}

	// Really only Halton samples from doubles, floats, or maybe shorts
	//    make sense.  Bools will likely look very unrandom, and I suspect
	//    bytes will only make sense for a very small number of samples.
	//
	// The first sample is sampIdx = 1.  (Not sampIdx = 0)
	double         dSample( unsigned int sampIdx );    /* in [0..1]     */
	float          fSample( unsigned int sampIdx );    /* in [0..1]     */
	unsigned short sSample( unsigned int sampIdx );    /* in [0..65535] */
	bool           bSample( unsigned int sampIdx );    /* true or false */
	unsigned char  cSample( unsigned int sampIdx );    /* in [0..255]   */

	// A pointer to a IGLUHalton1D could have type IGLUHalton1D::Ptr
	typedef IGLUHalton1D *Ptr;

private:
	// The base used for this set of Halton numbers
	unsigned int m_base;
	float        m_1_basef;  // A precomputed 1/base for float samples
	double       m_1_based;  // A precomputed 1/base for double samples
};


// End iglu namespace
}


#endif

