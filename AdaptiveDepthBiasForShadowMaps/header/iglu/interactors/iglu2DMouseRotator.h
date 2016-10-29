/************************************************************************/
/* iglu2DMouseRotator.h                                                 */
/* ---------------------                                                */
/*                                                                      */
/* This is a simple class that gives a rotation matrix around the       */
/*    z-axis, with the amount depending on the amount of rotation the   */
/*    mouse moves around a specified point on screen.                   */
/*                                                                      */
/* Inputs are:                                                          */
/*    * Integers width & height.  Often these define the window width & */
/*      height, they can also be thought of as the number of pixels you */
/*      need move the mouse to translate by xVec and yVec.              */
/*    * Integers x & y, specifying the point on screen to rotate around */
/*      these should be specified in WINDOW coordinates.                */
/*    * float speed, is a multipler controlling how fast mouse rotation */
/*      maps to world-space rotation.                                   */
/*                                                                      */
/* Chris Wyman (1/24/2012)                                              */
/************************************************************************/


#ifndef IGLU_2DMOUSE_ROTATOR_H
#define IGLU_2DMOUSE_ROTATOR_H

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

#include "iglu/interactors/igluMouseInteractor.h"

namespace iglu {

class vec2;

class IGLU2DMouseRotator : public IGLUMouseInteractor
{
public:

	IGLU2DMouseRotator( int width, int height, 
		                int xRotCenter, int yRotCenter, 
						float speed=1.0 );
	virtual ~IGLU2DMouseRotator()                                         {}

	// Functions for updating the interactor due to mouse input.           
	//    x & y are the window coordinates of the current mouse position. 
	//    The return value of UpdateOnMotion() is usually ignored, but should
	//    return something related how far we've moved (could be in pixels or
	//    degrees, for instance).
	virtual void SetOnClick( int x, int y );
	virtual float UpdateOnMotion( int x, int y );
	virtual void Release( void );

	// What to do when the interactor is reset
	virtual void Reset( void );

	// Update where the center of rotation is
	void ResetRotationCenter( int xPos, int yPos );         

	// A pointer to a IGLUMouseTranslator could have type IGLUMouseTranslator::Ptr
	typedef IGLU2DMouseRotator *Ptr;

private:
	int m_lastPos[2];     // Tracks last mouse position
	float m_speed;        // A multiplicative factor for how fast to move
	float m_curAngle;     // Total rotation seen so far
	vec3 m_rotCtr;        // The center point in window coords for rotation (z=0)
};


// End namespace iglu
}

#endif
