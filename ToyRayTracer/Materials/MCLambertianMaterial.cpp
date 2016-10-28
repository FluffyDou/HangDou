/******************************************************************/
/* MCLambertianMaterial.cpp                                       */
/* -----------------------                                        */
/*                                                                */
/* The file defines a material that uses for MCR                  */
/*                                                                */
/* Hang Dou (03/31/2012)                                          */
/******************************************************************/

#include "Materials/MCLambertianMaterial.h"
#include "Core/Scene.h"

MCLambertianMaterial::MCLambertianMaterial( const Color &matlColor ): matlColor(matlColor)
{
}


Color MCLambertianMaterial::Shade( Ray &ray, Scene &scn ) const
{
    // Double check that we have not bounced too many times
    if (ray.bounce > scn.bounceLimit)
    {
        return Color(0.0, 0.0, 0.0);
    }

    /********** record world position of certain intersection --- ray.bounce starts from 1 ( gets ++ in scene.Trace() ) *******/
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

    // Remember, we are doing Monte Carlo integration of the rendering equation:
    //    Int[ f_r(x,w_in,w_out) * L_in(x,w_in) * dot(w_in,pNorm) ]
    // For Lambertian BRDFs, f_r is a constant of albedo/pi, which is surfaceColor/pi.
    // Using Monte Carlo integration, the above integral becomes:
    //    (surfaceColor/(pi*N)) * Sum[ (L_in(x,w_in) * dot(w_in,pNorm)) / prob(w_in) ]
    // Here, N is the number of samples used.  In this shader's notation, w_in is called 
    // reflectVec.  reflectVec was selected uniformly over the hemisphere, so 
    // prob(w_in) = 1/(2*pi), which means this further simplifies to:
    //    (2*pi*surfaceColor/(pi*N)) * Sum[ L_in(x,reflectVec) * dot(reflectVec,pNorm) ]
    // Since we don't know N here, we let the camera do the divide by N.

    return matlColor * bounceResult * rndDir.Dot(ray.hitNormal) * 2.0;
}