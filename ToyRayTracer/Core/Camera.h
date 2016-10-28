/******************************************************************/
/* Camera.h                                                       */
/* -------------                                                  */
/*                                                                */
/* The file defines the camera class for the ray tracer.  The     */
/*     Important functions are the constructor, and the           */
/*     GenerateRay() method.  GenerateRay() takes pixel coords    */
/*     (either integer or pixel) and returns a vector from the    */
/*     eye (at the point 'eye') through that pixel in the image.  */
/*     Camera() takes in the eye point, look at point, an up      */
/*     vector, field of view, and the image resolution.           */
/*                                                                */
/* Chris Wyman (10/26/2006)                                       */
/******************************************************************/

#ifndef CAMERA_H
#define CAMERA_H 

#include "DataTypes/datatypes.h"
#include "Core/Ray.h"
#include "DataTypes/MathDefs.h"

#include <iostream>
#include <math.h>

class Scene;

class Camera {

  // Store the camera a way we like to think about it (and OpenGL does)
  vec3  eye, at, up;
  float fovy;

  // Store the screen size
  int screenWidth, screenHeight;

  // Store the camera in an easy-to-use format for ray tracing
  vec3 U, V, W;

public:

  // sets up a camera 
  Camera( const vec3 &eye, const vec3 &at, const vec3 &up, float fovy, 
	      int screenWidth, int screenHeight );

  // generate a ray based on previous set camera parameters 
  vec3 GenerateRay( float x, float y );

  // accessor functions 
  inline vec3 GetEye( void ) const { return eye; }
  inline int GetScreenWidth( void ) const { return screenWidth; }
  inline int GetScreenHeight( void ) const { return screenHeight; }

  void SetScreenRatio(int w, int h);
};


#endif

