/******************************************************************/
/* LightMaterial.h                                                */
/* -----------------------                                        */
/*                                                                */
/*                                                                */
/* Hang Dou (04/01/2012)                                          */
/******************************************************************/

#ifndef LIGHTMATERIAL_H
#define LIGHTMATERIAL_H

#include "DataTypes/datatypes.h"
#include "Materials/Material.h"

class LightMaterial : public Material {
private:
    Color matlColor;
public:
    LightMaterial( const Color &matlColor );
    virtual ~LightMaterial();
    virtual Color Shade( Ray &ray, Scene &scn ) const;
};

#endif
