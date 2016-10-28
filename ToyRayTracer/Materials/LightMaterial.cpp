/******************************************************************/
/* LightMaterial.cpp                                              */
/* -----------------------                                        */
/*                                                                */
/* The file defines a two bounce light material that uses for MCR */
/*                                                                */
/* Hang Dou (03/31/2012)                                          */
/******************************************************************/

#include "Materials/LightMaterial.h"
#include "Core/Scene.h"

LightMaterial::LightMaterial( const Color &matlColor ): matlColor(matlColor)
{
}


LightMaterial::~LightMaterial(void)
{
}

Color LightMaterial::Shade( Ray &ray, Scene &scn ) const
{
    // when hit the light, stop bouncing
    //ray.bounce = scn.bounceLimit+1; --- no need, because no bounce here

    //return intense = lightColor * dot( -ray.direction, lightNorm );
    //return Color( matlColor * MAX( 0.0, (ray.hitNormal).Dot(-1.0*ray.direction) ) );
    /********** record world position of certain intersection --- ray.bounce starts from 1 ( gets ++ in scene.Trace() ) *******/

    // Double check that we have not bounced too many times
    if (ray.bounce > scn.bounceLimit)
    {
        return matlColor;
    }

    (scn.X->worldSpaceCoord)[ray.bounce-1] = ray.GetHitPoint();

    /*********** record the surface normal of the certain intersection ****************/
    (scn.X->surfaceNorm)[ray.bounce-1] = ray.hitNormal;

    /*********** record the texture color of the certain intersection ****************/
    if(ray.bounce == 1)
        (scn.X->texValue)[0] = this->matlColor.ColorVec();

    // Offset the hit point slightly in the normal direction to help avoid self intersections
    vec3 offsetOrig = ray.GetHitPoint() + 0.01 * ray.hitNormal;
    // get the random direction over the hemisphere
    Random rnd( unsigned long( scn.passSampleSeed * abs(0.7*ray.direction.X()+0.3*ray.direction.Y()) ) );
    // global random seed ++
    scn.passSampleSeed += 1;

    float randomPara1 = 0.0;
    float randomPara2 = 0.0;
    vec3  rndDir = ( rnd.RandomHemisphereVector(ray.hitNormal, randomPara1, randomPara2) ).vNormalize();

    /************ record parameter for bouncing direction  of the certain intersection *********/
    (scn.X->ranPara)[ray.bounce-1] = vec2(randomPara1, randomPara2);

    // bounce the ray into a random direction
    Ray bounceRay( offsetOrig, rndDir, ray.bounce );
    Color bounceResult = scn.TraceRay(bounceRay); // after trace the ray, ray.bounce will increment

    return matlColor +  bounceResult * rndDir.Dot(ray.hitNormal) * 2.0;
    //return matlColor;
}