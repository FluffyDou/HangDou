/************************************************************************/
/* igluCameraRotator.h                                                  */
/* ------------------                                                   */
/*                                                                      */
/* This class implements a simple rotator around the eye position, used */
/*    for mouse rotations of the camera view point.                     */
/*                                                                      */
/* Inputs are:                                                          */
/*    * Integers width & height.  To rotate in a plausible manner, the  */
/*      rotator needs the correct window (or subwindow) width & height. */
/*                                                                      */
/* Chris Wyman (01/14/2013)                                             */
/************************************************************************/


#ifndef IGLU_CAMERA_ROTATOR_H
#define IGLU_CAMERA_ROTATOR_H

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

#include "iglu/interactors/igluMouseInteractor.h"

namespace iglu {

class IGLUCameraRotator : public IGLUMouseInteractor
{
public:

	/* sets up a trackball with window size of width x height */
	IGLUCameraRotator( int width, int height );
	virtual ~IGLUCameraRotator()                                  {}

	// Functions for updating the interactor due to mouse input.           
	//    x & y are the window coordinates of the current mouse position. 
	virtual void SetOnClick( int x, int y );
	virtual float UpdateOnMotion( int x, int y );
	virtual void Release( void );

	// What to do when the interactor is reset
	virtual void Reset( void );

	// A pointer to a IGLUCameraRotator could have type IGLUCameraRotator::Ptr
	typedef IGLUCameraRotator *Ptr;

private:
	float m_lastPos[3]; 
};


// End namespace iglu
}

#endif
