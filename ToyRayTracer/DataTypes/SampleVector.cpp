/******************************************************************/
/* SampleVector.cpp                                               */
/* -----------------------                                        */
/*                                                                */
// This pool contains all sample vectors                          */
/*                                                                */
/* Hang Dou (04/12/2012)                                          */
/******************************************************************/

#include "Datatypes/SampleVector.h"


SampleVector::SampleVector()
{
    screenPos[0]       = vec2(0.0, 0.0);
    sampleColor[0]     = vec3(0.0, 0.0, 0.0);
    texValue[0]        = vec3(0.0, 0.0, 0.0);
    ranPara[0]         = vec2(0.0, 0.0);
    ranPara[1]         = vec2(0.0, 0.0);
    worldSpaceCoord[0] = vec3(0.0, 0.0, 0.0);
    worldSpaceCoord[1] = vec3(0.0, 0.0, 0.0);
    surfaceNorm[0]     = vec3(0.0, 0.0, 0.0);
    surfaceNorm[1]     = vec3(0.0, 0.0, 0.0);

}


SampleVector::~SampleVector()
{  
}


void SampleVector::Format()
{
    screenPos[0]       = vec2(0.0, 0.0);
    sampleColor[0]     = vec3(0.0, 0.0, 0.0);
    texValue[0]        = vec3(0.0, 0.0, 0.0);
    ranPara[0]         = vec2(0.0, 0.0);
    ranPara[1]         = vec2(0.0, 0.0);
    worldSpaceCoord[0] = vec3(0.0, 0.0, 0.0);
    worldSpaceCoord[1] = vec3(0.0, 0.0, 0.0);
    surfaceNorm[0]     = vec3(0.0, 0.0, 0.0);
    surfaceNorm[1]     = vec3(0.0, 0.0, 0.0);

}