/******************************************************************/
/* SampleVector.h                                                 */
/* -----------------------                                        */
/*                                                                */
// This pool contains all sample vectors                          */
/*                                                                */
/* Hang Dou (04/12/2012)                                          */
/******************************************************************/

#ifndef SAMPLEVECTOR_H
#define SAMPLEVECTOR_H

#include "DataTypes/vec2.h"
#include "DataTypes/vec3.h"
#include "DataTypes/MathDefs.h"
#include "DataTypes/Color.h"

class SampleVector
{
public:
    vec2  screenPos[1];
    vec3  sampleColor[1];
    vec3  texValue[1];
    vec2  ranPara[2];
    vec3  worldSpaceCoord[2];
    vec3  surfaceNorm[2];
    
public:
    SampleVector();

    ~SampleVector();

    void Format();
};



#endif