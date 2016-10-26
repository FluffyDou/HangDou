/*********************************************************************/
// CPUTimer.cpp
// A high resolution (i.e., nanosecond) timer class on windows. This uses the CPU's internal clock-cycle count
// to compute ultra-high resolution timings.  Note that times may not be updated every nanosecond, depending on
// how often the OS queries the CPU for the cycle count.
// This class is basically grabbed from IGLU

// Hang 10/08/2014
/*********************************************************************/

#include "CPUTimer.h"

#pragma comment(lib, "kernel32.lib")
//typedef LARGE_INTEGER TimerStruct;

void CPUTimer::GetHighResolutionTime( LARGE_INTEGER *t ) 
{
    QueryPerformanceCounter( t ); 
}

double CPUTimer::ConvertTimeDifferenceToSec( LARGE_INTEGER *end, LARGE_INTEGER *begin ) 
{
    LARGE_INTEGER freq;
    QueryPerformanceFrequency( &freq ); 
    return (end->QuadPart - begin->QuadPart)/(double)freq.QuadPart; 
}

double CPUTimer::ConvertTimeDifferenceToMSec( LARGE_INTEGER *end, LARGE_INTEGER *begin ) 
{
    LARGE_INTEGER freq;
    QueryPerformanceFrequency( &freq );
    return (end->QuadPart - begin->QuadPart)/(1e-3*(double)freq.QuadPart); 
}

CPUTimer::CPUTimer(void)
{
    GetHighResolutionTime( &lastTime );
}

void CPUTimer::Start( void )
{
    GetHighResolutionTime( &lastTime ); 
}

double CPUTimer::GetTime( void )
{
    GetHighResolutionTime( &curTime );  
    return ConvertTimeDifferenceToMSec( &curTime, &lastTime ); 
}


double CPUTimer::Tick( void )
{
    curTime = lastTime; 
    GetHighResolutionTime( &lastTime ); 
    return ConvertTimeDifferenceToMSec( &lastTime, &curTime ); 
}

CPUTimer::~CPUTimer(void)
{
}
