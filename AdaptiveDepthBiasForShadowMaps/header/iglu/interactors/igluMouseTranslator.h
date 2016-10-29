/************************************************************************/
/* igluMouseTranslator.h                                                */
/* ---------------------                                                */
/*                                                                      */
/* This is a simple class that implements an translator based on mouse  */
/*    interactions (similar to the IGLUTrackball class).                */
/*                                                                      */
/* Inputs are:                                                          */
/*    * Integers width & height.  Often these define the window width & */
/*      height, they can also be thought of as the number of pixels you */
/*      need move the mouse to translate by xVec and yVec.              */
/*    * vec3 xVec.  The direction the matrix should translate when the  */
/*      user moves the mouse towards the right of the screen.           */
/*    * vec3 yVec.  The direction the matrix should translate when the  */
/*      user moves the mouse towards the top of the screen.             */
/*    * float speed.  A multiplier.  By default, moving <width> pixels  */
/*      translates by xVec.  With speed=2, it would translate by 2*xVec */
/*                                                                      */
/* Defaults are given for xVec, yVec, and speed, though depending on    */
/* the usage, these may cause either huge or unnoticable translations.  */
/*                                                                      */
/* Chris Wyman (1/24/2012)                                              */
/************************************************************************/


#ifndef IGLU_MOUSE_TRANSLATOR_H
#define IGLU_MOUSE_TRANSLATOR_H

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

#include "iglu/interactors/igluMouseInteractor.h"

namespace iglu {

class IGLUMouseTranslator : public IGLUMouseInteractor
{
public:

	// Sets up an interactor the window size, a multiplicative speed,
	//    and the directions the matrix should translate in for movements
	//    in the x-axis and the y-axis.
	// Movements of <width>  pixels in X cause a translation of speed*xVec
	// Movements of <height> pixels in Y cause a translation of speed*yVec
	IGLUMouseTranslator( int width, int height, float speed=1.0, 
		                 const vec3& xVec=vec3::XAxis(), 
						 const vec3& yVec=vec3::YAxis() );
	virtual ~IGLUMouseTranslator()                                         {}

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

	// A pointer to a IGLUMouseTranslator could have type IGLUMouseTranslator::Ptr
	typedef IGLUMouseTranslator *Ptr;

private:
	int m_lastPos[2];     // Tracks last mouse position
	float m_speed;        // A multiplicative factor for how fast to move
	vec3 m_xVec, m_yVec;  // Directions to translate for mouse moves in x & y
};


// End namespace iglu
}

#endif
