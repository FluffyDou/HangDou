/******************************************************************/
/* AmbientOcclusionMaterial.h                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines a simple ambient occlusion material using     */
/*    random sampling.                                            */
/*                                                                */
/* Chris Wyman (08/30/2010)                                       */
/******************************************************************/

#ifndef AMBIENTOCCLUSIONMATERIAL_H
#define AMBIENTOCCLUSIONMATERIAL_H

#include "DataTypes/datatypes.h"
#include "Materials/Material.h"

class AmbientOcclusionMaterial : public Material {
private:
	int sqrtNumSamples;
public:
	AmbientOcclusionMaterial( int sqrtNumSamples );
	virtual Color Shade( Ray &ray, Scene &scene ) const;
};

#endif

