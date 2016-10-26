/*********************************************************************/
// CPUTimer.h
// A high resolution (i.e., nanosecond) timer class on windows. This uses the CPU's internal clock-cycle count
// to compute ultra-high resolution timings.  Note that times may not be updated every nanosecond, depending on
// how often the OS queries the CPU for the cycle count.
// This class is basically grabbed from IGLU

// Hang 10/08/2014
/*********************************************************************/

#pragma once

#include <windows.h>

class CPUTimer
{
public:

    CPUTimer(void);

    ~CPUTimer(void);

    // Start / restart timing
    void Start( void );

    // Get the current time (in ms), without resetting the timer
    double GetTime( void );     

    // Get the current time (in ms), and restart the timer
    double Tick( void );

    // A pointer to a CPUTimer could have type CPUTimer::Ptr
    typedef CPUTimer *Ptr;

private:
    LARGE_INTEGER lastTime, curTime;

private:
    void   GetHighResolutionTime( LARGE_INTEGER *t );
    double ConvertTimeDifferenceToSec( LARGE_INTEGER *end, LARGE_INTEGER *begin );
    double ConvertTimeDifferenceToMSec( LARGE_INTEGER *end, LARGE_INTEGER *begin );

};

