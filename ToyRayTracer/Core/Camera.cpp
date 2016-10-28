/******************************************************************/
/* Camera.cpp                                                     */
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

#include "Camera.h"
#include "Scene.h"
#include "DataTypes/MathDefs.h"

Camera::Camera( const vec3 &eye, const vec3 &at, const vec3 &up, float fovy, 
	  int screenWidth, int screenHeight ) : eye(eye), at(at), up(up), fovy(fovy),
	  screenWidth(screenWidth), screenHeight(screenHeight)
{
    W = vec3(at - eye).vNormalize();
    U = W.Cross( up ).vNormalize();
    V = U.Cross( W ).vNormalize();

    //float vLen = tanf( DEG2RAD( fovy * 0.5 ) );
    float vLen = 0.0175; // height is 0.025
    float ulen = vLen * ( (float)screenWidth / (float)screenHeight );

    U = U * ulen;
    V = V * vLen;
}


// generate a ray direction for pixel x,y
vec3 Camera::GenerateRay( float x, float y )
{

    vec2 d = vec2((float)x, (float)y)/vec2((float)screenWidth, (float)screenHeight) * 2.0f - vec2(1.0, 1.0);

    return vec3(U * d.X() + V * d.Y() + 0.035*W).vNormalize(); // the focal length is 0.035
}

void Camera::SetScreenRatio(int w, int h)
{
    screenWidth  = w;
    screenHeight = h;
}

