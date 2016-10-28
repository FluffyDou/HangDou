/******************************************************************/
/* AmbientOcclusionMaterial.cpp                                   */
/* -----------------------                                        */
/*                                                                */
/* The file defines a material class that implements a very       */
/*     simple ambient occlusion (a darkening when objects are     */
/*     near to each other).                                       */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "Materials/AmbientOcclusionMaterial.h"
#include "Core/Ray.h"
#include "Objects/Primitive.h"
#include "Core/Scene.h"

AmbientOcclusionMaterial::AmbientOcclusionMaterial( int sqrtNumSamples ) :
	sqrtNumSamples(sqrtNumSamples)
{
}


Color AmbientOcclusionMaterial::Shade( Ray &ray, Scene &scene ) const
{
	// Get the ray's hit position
	vec3 hitPos = ray.GetHitPoint();
	
	// We'll send off ambient occlusion feelers in numerous random directions
	//   and count how many hit other objects
	int visibilityCount = 0;

	// Actually send the rays
	for (int i=0; i < sqrtNumSamples*sqrtNumSamples; i++)
	{
		// Compute a random direction for the ray
		vec3 randomDir = scene.rng->RandomSphereVector();
		
		// Setup our ray and shoot it
		Ray rndRay( hitPos, randomDir );
		if ( !scene.VisibilityRay( rndRay ) )
			visibilityCount++;
	}

	// The count of visibile rays is not the actual color.  We need to take a
	//   percentage.  Then becuase I implemented this in a (stupidly) simple way
	//   I also have to multiply by 2.
	float grayCol = (float)visibilityCount/(0.5f*sqrtNumSamples*sqrtNumSamples);

	// Return our final color!
	return Color( grayCol, grayCol, grayCol );
}

