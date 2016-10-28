/******************************************************************/
/* raytracer.h                                                    */
/* -----------------------                                        */
/*                                                                */
/* The basic include file for the main raytacing files.  This     */
/*     includes most of the relevent material and object headers  */
/*     and defines prototypes for useful functions.               */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <iostream>
#include <stdio.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <time.h>

#include "Utils/Random.h"
#include "Utils/frameRate.h"
#include "DataTypes/datatypes.h"
#include "Core/Camera.h"
#include "Core/Scene.h"
#include "Objects/Sphere.h"
#include "Objects/InfinitePlane.h"
#include "Objects/Triangle.h"
#include "Objects/Cylinder.h"
#include "Materials/ConstantColorMaterial.h"
#include "Materials/LambertianMaterial.h"
#include "Materials/MCLambertianMaterial.h"
#include "Materials/AmbientOcclusionMaterial.h"
#include "Materials/LightMaterial.h"
#include "DataTypes/MathDefs.h"

#include "Shaders/texShader.h"
//#include "GLObjects/frameBuf.h"

#include "iglu.h"

#include "glm.h"

#ifndef myDebug
#define myDebug(var) std::cout<<var<<"\n";
#endif

const int LINEAR          = 0;
const int NEAREST         = 1;

const int SCREENPOS       = 1;
const int RANPARA         = 2;
const int WORLDSPACECOORD = 3;
const int SURFACENORMAL   = 4;
const int TEXVALUE        = 5;
const int SAMPLECOLOR     = 6;
const int RENDERRESULT    = 7;

const int SCREENWIDTH  = 600;
const int SCREENHEIGHT = 600;

const float FOV       =  60.0;
const vec3  LOOKEYE   =  vec3( 278.0, 273.0, -500.0);
const vec3  LOOKAT    =  vec3( 278.0, 273.0, -499.0);
const vec3  LOOKUP    =  vec3(0.0, 1.0, 0.0);

const int SPP         = 16; // sample limit per pixel

int sampleCount        = 0; // sample pass count
int stageFlag1         = 0; // post processing flag --- avoid processing the data every frame
int renderFlag         = 0; // choose to render which vector element of the sample vector

// a scene contains the geometry and the camera
Scene* myScene = 0;
// a CPU memory to store ray tracing value
Image* myImage = 0;
// a CPU memory to store the average value of sample vectors
Image* featureImage = 0;

// frame buffer object
//frameBuf*  renderBuffer  = 0;

// shader to display data in a texture
texShader* textureShader = 0;
// A frame rate counter
FrameRate *fRate = 0;

void readInObjFile(Group* grp, Material** mt);
void createCornellBox(Group* grp, Material** mt);                                  // create a cornell box
Scene *SceneSetup( void );                                                         // set up a scene including geometry, material and camera
void ShootRays( Image &image, Scene *scene, int samplePass );                      // trace the rays and store the data into a texture
//frameBuf* RenderBufInit(int screenWidth, int screenHeigth);                      // initialize the Render Texture for frame buffer object
texShader* TexShaderSetup(char* vshader, char* fshader);                           // initialize Shader

void AverageSampleVectors( Image &destination, Image &source );                     // average the sample vector values to render them out
void UpdateCanvus( Image &image, int vectorFlag );                                                  // update the texture of Image for rendering 

void myIdle();                                  // GL idle function
void myDisplay();                               // GL display function
void myReshape(int w, int h);                   // reshape function
void myKeys(unsigned char key, int x, int y);   // GL key call back function

// Displays a timer on the lower left corner of the screen
void DisplayTimer( float fps );

// This disables a particularly stupid, Microsoft-only compiler warning.
#pragma warning( disable:4996 )

#endif

