/******************************************************************/
/* Scene.h                                                        */
/* -----------------------                                        */
/*                                                                */
/* The file defines a scene class that encapsulates all the       */
/*     information necessary to render an image with a ray tracer */
/* Also note that this class includes the TraceRay() method,      */
/*     which actually traces a ray through the scene.             */
/* Scenes should be setup similar to the example(s) in the file   */
/*     SceneSetup.cpp.                                            */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef SCENE_H
#define SCENE_H

#include "DataTypes/Color.h"
#include "Objects/Group.h"
#include "Objects/Object.h"
#include "DataTypes/MathDefs.h"
#include "Core/Camera.h"
#include "DataTypes/Array1D.h"
#include "Utils/Random.h"
#include "Lights/Light.h"
#include "Datatypes/SampleVector.h"

class Ray;

class Scene {
public:
	Camera *camera;           // Scene camera.  There is only one (for now).
	Object *geometry;         // Scene geometry.  There is only one, as it should be a container object.
    Light  *light;            // one light source
	Color   backgroundColor;  // A background color. 
	Random *rng;			  // A random number generator
    SampleVector *X;          // a sample vector containing all the features we get from a ray tracing

    vec3 sceneAmbient;    //= vec3(0.3, 0.3, 0.3);
    vec3 materialAmbient; // = vec3(0.3, 0.3, 0.3);

    int bounceLimit;
    int passSampleSeed;

public:
	// Constructors & destructors
	Scene();						 
	~Scene();

	// Get the width and height of the image
	inline int GetWidth( void ) { return camera->GetScreenWidth(); }
	inline int GetHeight( void ) { return camera->GetScreenHeight(); }
	
	// Trace a ray through the scene and return a color
	Color TraceRay( Ray &r );

	// Trace a ray through the scene and see if it hit anything
	bool VisibilityRay( Ray &r );

};



#endif

