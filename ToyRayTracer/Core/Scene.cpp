/******************************************************************/
/* Scene.cpp                                                      */
/* -----------------------                                        */
/*                                                                */
/* The file defines a scene class that encapsulates all the       */
/*     information necessary to render an image with a ray tracer */
/* Also not that this class includes the TraceRay() method, which */
/*     actually traces a ray through the scene.                   */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#include "DataTypes/Color.h"
#include "Core/Scene.h"
#include "Core/Camera.h"
#include "Objects/Group.h"
#include "Objects/Primitive.h"
#include "Materials/Material.h"
#include "Materials/LambertianMaterial.h"
#include "Materials/ConstantColorMaterial.h"
#include "Objects/Sphere.h"
#include "Lights/Light.h"

Scene::Scene()
{
	rng = new Random();
    passSampleSeed = 0;
}

// Free all the allocated memory inside the scene.
Scene::~Scene()
{
	if (camera)     delete camera;
	if (geometry)   delete geometry;
	if (rng)        delete rng;
    if (X)          delete X;
}

// Intersect a ray with the scene, return it's color
Color Scene::TraceRay( Ray &r )
{
	// Find the closest intersection in the geometry
	geometry->Intersect( r );

	// If there we no intersection, return the background color
	if (!r.WasIntersection())
		return backgroundColor;

    // Else, we hit something
    r.bounce++;

	// Get the shade using the object's material type
    return r.hitObj->GetMaterial()->Shade( r, *this );
}


// Intersect a ray with the scene, return true if an intersection occurred
bool Scene::VisibilityRay( Ray &r )
{
	// Find the closest intersection in the geometry
	geometry->Intersect( r );

	// Return if there we an intersection
	return r.WasIntersection();
}



