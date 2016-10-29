/************************************************************************/
/* igluTrackball.h                                                      */
/* ------------------                                                   */
/*                                                                      */
/* This is a class that implements a simple virtual trackball for       */
/*    mouse-based 3D rotations.  The trackball is not a very good       */
/*    interface, but it is quick and easy to use and allows for simple  */
/*    addition of UI-based rotations.                                   */
/*                                                                      */
/* I have been using this trackball for over 10 years now, and this     */
/*    class is nearly identical to the Trackball.cpp class I used for   */
/*    most of those years.  This one is different in that it derives    */
/*    from a base mouse-interactor class to allow sharing of interfaces */
/*    between other simple mouse interactions.  The only difference     */
/*    in the class is that method names including 'Trackball' have been */
/*    renamed suitable for a more generic interface.                    */
/*                                                                      */
/* Inputs are:                                                          */
/*    * Integers width & height.  To rotate in a plausible, 'correct'   */
/*      manner, the trackball needs the window width & height.          */
/*                                                                      */
/* Chris Wyman (01/24/2012)                                             */
/************************************************************************/


#ifndef IGLU_TRACKBALL_H
#define IGLU_TRACKBALL_H

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>

#include "iglu/interactors/igluMouseInteractor.h"

namespace iglu {

class IGLUTrackball : public IGLUMouseInteractor
{
public:

	/* sets up a trackball with window size of width x height */
	IGLUTrackball( int width, int height );
	virtual ~IGLUTrackball()                                              {}

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

	// A pointer to a IGLUTrackball could have type IGLUTrackball::Ptr
	typedef IGLUTrackball *Ptr;

private:
	float m_lastPos[3]; 

	void trackball_ptov(int x, int y, int width, int height, float v[3]);
};


// End namespace iglu
}

#endif
