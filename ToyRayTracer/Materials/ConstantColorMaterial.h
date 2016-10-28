/******************************************************************/
/* ConstantColorMaterial.h                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a basic material class that simply colors     */
/*     the geometry a constant color.                             */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef CONSTANTCOLORMATERIAL_H
#define CONSTANTCOLORMATERIAL_H

#include "DataTypes/datatypes.h"
#include "Materials/Material.h"

class ConstantColorMaterial : public Material {
private:
	Color matlColor;
public:
	ConstantColorMaterial( const Color &matlColor );
	virtual ~ConstantColorMaterial() {}
	virtual Color Shade( Ray &ray, Scene &scene ) const;
};

#endif

