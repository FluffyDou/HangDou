/******************************************************************/
/* LembertianMaterial.h                                           */
/* -----------------------                                        */
/*                                                                */
/* The file defines a material that uses the angle between the    */
/*     light and the surface normal to compute a simple smooth    */
/*     shaded color along the surface.                            */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef LAMBERTIANMATERIAL_H
#define LAMBERTIANMATERIAL_H

#include "DataTypes/datatypes.h"
#include "Materials/Material.h"

class LambertianMaterial : public Material {
private:
	Color matlColor;
public:
	LambertianMaterial( const Color &matlColor );
	virtual ~LambertianMaterial() {}
	virtual Color Shade( Ray &ray, Scene &scn ) const;
};

#endif

