/******************************************************************/
/* Material.h                                                     */
/* -----------------------                                        */
/*                                                                */
/* The file defines a base material class that is inherited by    */
/*     more complex material types.                               */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef MATERIAL_H
#define MATERIAL_H

#include "DataTypes/datatypes.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma warning( disable:4996 )

class Ray;
class Scene;

class Material {
public:
	Material() {}
	virtual ~Material() {}
	virtual Color Shade( Ray &ray, Scene &scene ) const = 0;
};


#endif

