/******************************************************************/
/* ConstantColorMaterial.cpp                                      */
/* -------------------------                                      */
/*                                                                */
/* The file defines a basic material class that simply colors     */
/*     the geometry a constant color.                             */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/


#include "Materials/ConstantColorMaterial.h"
#include "Core/Scene.h"
#include "DataTypes/MathDefs.h"

ConstantColorMaterial::ConstantColorMaterial( const Color &matlColor ) :
	matlColor(matlColor)
{
}

Color ConstantColorMaterial::Shade( Ray &ray, Scene &scene ) const
{
	// Shading constant color materials is tough, let me tell you!
	return matlColor;
}

