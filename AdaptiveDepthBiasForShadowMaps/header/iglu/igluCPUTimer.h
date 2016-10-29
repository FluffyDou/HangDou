/************************************************************************/
/* igluCPUTimer.h                                                       */
/* --------------                                                       */
/*                                                                      */
/* A high resolution (i.e., nanosecond) timer class based on system-    */
/* dependent OS calls.  This uses the CPU's internal clock-cycle count  */
/* to compute ultra-high resolution timings.  Note that times may not   */
/* be updated every nanosecond, depending on how often the OS queries   */
/* the CPU for the cycle count.                                         */
/*                                                                      */
/* System specific info:                                                */
/*     * Windows:  Should work out-of-the-box by including this header  */
/*                 and using the IGLUTimer class                        */
/*     * Linux:  Must link with the "rt" library to use this class      */
/*                 (i.e., add "-lrt" to the command line when linking)  */
/*     * MacOS:  Must link with the "CoreServices" framework to use the */
/*                 class (add "-framework CoreServices" to the g++      */
/*                 command line while linking)                          */
/*                                                                      */
/* Chris Wyman (06/23/2011)                                             */
/************************************************************************/

#ifndef __IGLU_TIMER_H__
#define __IGLU_TIMER_H__

// The CPU timer needs an appropriate OS-dependent timer structure.  
//    Set that up here.  Other OS-dependent stuff is handled in the library.

/************************************************************************/
/*  The CPU timer needs an appropriate OS-dependent structure to store  */
/*      intermediate timing information.  Set that up here.  Other OS-  */
/*      specific stuff is actually handled in the library.              */
/************************************************************************/

#if defined(WIN32) && defined(_MSC_VER) && !defined(USING_MSVC)
	#include <windows.h>
namespace iglu {
    typedef LARGE_INTEGER TimerStruct;
}

#elif defined(__APPLE__) && !defined(USING_MACOSX) && !defined(USING_LINUX) 
	// The following #includes occur in the library.  
    //      One may be needed for the uint64_t(?)  We do not test regularly 
    //      on MacOS.
	//#include <CoreServices/CoreServices.h>
    //#include <mach/mach.h>
    //#include <mach/mach_time.h>
    
namespace iglu {
    typedef uint64_t TimerStruct;
}

#elif defined(__GNUC__) && !defined(USING_MACOSX) && !defined(USING_LINUX) 
	#include <time.h>
namespace iglu {
    typedef struct timespec TimerStruct;
}

#endif


/************************************************************************/
/*  Now define our timer class.  This is the important part!            */
/************************************************************************/

namespace iglu {

class IGLUCPUTimer
{
private:
	TimerStruct lastTime, curTime;
public:
	// Allocate / start a timer
	IGLUCPUTimer();

	// Start / restart timing
	void Start( void );

	// Get the current time (in ms), without resetting the timer
	double GetTime( void );     

	// Get the current time (in ms), and restart the timer
	double Tick( void );

	// Get the current time (in ms), and "stop" the timer
	//    (This is the same as IGLUCPUTimer::Tick(), but unifies the interface
	//     between IGLUCPUTimer and IGLUGPUTimer [where Tick ()and Stop() vary])
	double Stop( void );

	// A pointer to a IGLUCPUTimer could have type IGLUCPUTimer::Ptr
	typedef IGLUCPUTimer *Ptr;
};

}

#endif


