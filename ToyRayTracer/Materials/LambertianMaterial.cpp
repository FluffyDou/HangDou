/******************************************************************/
/* LembertianMaterial.cpp                                         */
/* -----------------------                                        */
/*                                                                */
/* The file defines a material that uses the angle between the    */
/*     light and the surface normal to compute a simple smooth    */
/*     shaded color along the surface.                            */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/


#include "Materials/LambertianMaterial.h"
#include "Core/Scene.h"
//#include "DataTypes/MathDefs.h"

LambertianMaterial::LambertianMaterial( const Color &matlColor ) :
	matlColor(matlColor)
{
}

Color LambertianMaterial::Shade( Ray &ray, Scene &scn ) const
{
    // Get the ray's hit position
    vec3 hitPos     = ray.GetHitPoint();
    float radiusSqr = 1.0;

    vec3 lightDir = (scn.light->GetPosition() - hitPos).vNormalize();

    vec3 ambientLighting = scn.sceneAmbient *matlColor.ColorVec();
    // 1/pi is the BRDF for lambertian
    vec3 diffuseLighting = ( scn.light->GetIntense() * ( 1.0/(4.0*M_PI*radiusSqr) ) * matlColor.ColorVec() ) * MAX( 0.0, (ray.hitNormal).Dot(lightDir) ) * (1.0/M_PI);
    //( 1.0/(4.0*M_PI*radiusSqr) ) *
    return Color( ambientLighting + diffuseLighting );
}

