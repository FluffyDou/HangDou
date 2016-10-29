/******************************************************************/
/* igluVector.h                                                   */
/* -----------------------                                        */
/*                                                                */
/* The file includes the various GPU vector data types that are   */
/*    included in IGLU.  These can be used for GLSL/OpenGL        */
/*    interfacing or as points/vectors in 1D-3D spaces directly.  */
/*                                                                */
/* Eventually, all these vector classes will behave as primitive  */
/*    datatypes with all appropriate mathematical operations and  */
/*    casting operations built in.  However, this functionality   */
/*    is broadly implemented as needed, so not all operations may */
/*    be currently overloaded.  Beware.                           */
/*                                                                */
/* Chris Wyman (09/28/2011)                                       */
/******************************************************************/

#ifndef IGLU_VECTOR_DATATYPES_H
#define IGLU_VECTOR_DATATYPES_H

// Vectors/scalers using half-floats
#include "vectors/half.h"


// Vectors using integers
#include "vectors/int2.h"
#include "vectors/int3.h"
#include "vectors/int4.h"

// Vectors using unsigned integers
#include "vectors/uint2.h"
#include "vectors/uint3.h"
#include "vectors/uint4.h"

// Vectors using floats
#include "vectors/vec2.h"
#include "vectors/vec3.h"
#include "vectors/vec4.h"

// Vectors using doubles
#include "vectors/dvec2.h"
#include "vectors/dvec3.h"
#include "vectors/dvec4.h"



#endif
