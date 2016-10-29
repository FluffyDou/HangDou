/******************************************************************/
/* igluSampling.h                                                 */
/* -----------------------                                        */
/*                                                                */
/* The file contains a list of sampling routines in IGLU, and     */
/*    explicitly includes the appropriate class headers.          */
/*                                                                */
/* Chris Wyman (09/30/2011)                                       */
/******************************************************************/

#ifndef IGLUSAMPLING_H
#define IGLUSAMPLING_H

// For Quasi-random sampling using the Halton sequence.
#include "sampling/igluHalton1D.h"

// Base class for 2D samplers
#include "sampling/igluSampler2D.h"

// Some 2D samplers
#include "sampling/igluHalton2D.h"
#include "sampling/igluHammersley2D.h"

// Indexing schemes
#include "sampling/igluSpaceFillingCurve2D.h"
#include "sampling/igluMortonCurve2D.h"

#endif

