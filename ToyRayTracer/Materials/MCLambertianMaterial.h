/******************************************************************/
/* MCLambertianMaterial.cpp                                       */
/* -----------------------                                        */
/*                                                                */
/* The file defines a material that uses for MCR                  */
/*                                                                */
/* Hang Dou (03/31/2012)                                          */
/******************************************************************/

#ifndef MCLAMBARTIANMATERIAL_H
#define MCLAMBARTIANMATERIAL_H

#include "DataTypes/datatypes.h"
#include "Materials/Material.h"
#include "Utils/Random.h"

class MCLambertianMaterial : public Material
{
    Color matlColor;

public:
    MCLambertianMaterial( const Color &matlColor );
    virtual Color Shade( Ray &ray, Scene &scn ) const;

};

#endif

