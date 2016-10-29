/************************************************************************/
/* igluGPUTimer.h                                                       */
/* --------------                                                       */
/*                                                                      */
/* A basic timer for determining GPU usage based on GPU cycle-counting  */
/*     hardware and the OpenGL timer-query extension/API.               */
/*                                                                      */
/* Chris Wyman (9/26/2011)                                              */
/************************************************************************/


#ifndef IGLUGPUTimer_H__
#define IGLUGPUTimer_H__

#include <GL/glew.h>

namespace iglu {

class IGLUGPUTimer
{
private:
	GLuint64 nsElapsed;
	GLuint64 totalElapsed;
	GLuint   timeQuery;
	bool     initialized;
public:
	// Allocate / start a timer
	IGLUGPUTimer();					

	// Start / restart timing
	void Start( void );

	// Get the current time (in ms), without resetting the timer
	double GetTime( void );

	// Get the current time (in ms), and restart the timer
	double Tick( void );

	// Get the current time (in ms), and stop the timer
	double Stop( void );

	// A pointer to a IGLUGPUTimer could have type IGLUGPUTimer::Ptr
	typedef IGLUGPUTimer *Ptr;
};


}

#endif


