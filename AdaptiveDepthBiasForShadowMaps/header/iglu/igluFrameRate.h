/************************************************************************/
/* igluFrameRate.h                                                      */
/* ------------------                                                   */
/*                                                                      */
/* This is a pretty simple class that allows you to start timing at the */
/*    beginning of a frame's rendering and display the current          */
/*    framerate at the end of the frame.  Additionally, because 1-frame */
/*    counters are notoriously noisy (due to context switches, etc),    */
/*    this class allows the averaging over multiple frames, specified   */
/*    to the constructor.  Please note this class works only if the     */
/*    screen continuously updates (e.g., the idle function calls        */
/*    glutPostRedisplay()) -- after all a "framerate" makes no sense    */
/*    otherwise.                                                        */
/*                                                                      */
/* Chris Wyman (12/7/2007)                                              */
/************************************************************************/


#ifndef IGLU__FRAMERATE_H__
#define IGLU__FRAMERATE_H__

namespace iglu {

class IGLUFrameRate
{
public:
	// Initializes the frame rate counter.  The input will be the number
	//    of previous frames that framerate will be averaged over.
	IGLUFrameRate( int avgOverFrames = 5 );

	// Free allocated memory
	~IGLUFrameRate() { if (timeArray) free( timeArray ); }

	// Call at the very beginning of a frame
	void StartFrame( void );

	// Call at the end of a frame.  The value returned is the framerate
	//    (i.e., frames per second) averaged over the last N frames, where 
	//    N is the value specified in the constructor.  
	float EndFrame( void );

	// This returns the same value as the last call to EndFrame().  If you
	//    have not yet called EndFrame(), the return value is 0.
	inline float GetLastFrameRate( void )  { return lastFrameRate; }

	// A pointer to a IGLUFrameRate could have type IGLUFrameRate::Ptr
	typedef IGLUFrameRate *Ptr;

private:
	int avgOverFrames;

	int currFrameCount;
	int frameIdx;
	float lastFrameRate;   // the last value computed in EndFrame().  Set to 0 initially.
	
	TimerStruct *timeArray;
};


// End iglu namespace
}

#endif


