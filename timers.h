#ifndef _TIMERS_INC_
#define _TIMERS_INC_

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#include <string.h>
#endif

class Timers
{
	// cumulative timer class for profiling
	int numTimers;
	double *timers;

#ifdef _WIN32
	double invTimerFreqMicrosec;
#endif
	double microsecs( void )
	{
#ifdef _WIN32
		LARGE_INTEGER timeValue;
		QueryPerformanceCounter(&timeValue);
		return timeValue.QuadPart * invTimerFreqMicrosec;
#else
		timeval t;
		gettimeofday( &t, 0 );
		return (double)(t.tv_sec*1e6 + t.tv_usec);
#endif
	}

public:
	Timers( int nTimers )
	{
		numTimers = nTimers;
		timers = new double[numTimers];
		memset(timers, 0, numTimers*sizeof(double));
#ifdef _WIN32
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		invTimerFreqMicrosec = 1.0e6 / (double)freq.QuadPart;
#endif

	}
	void Clear( void )
	{
		for ( int i = 0; i < numTimers; i++ )
			timers[i] = 0.0;
	}

	void ClearTimer(int timer)
	{
		timers[timer] = 0.0;
	}

	void StartTimer( int timer )
	{
		timers[timer] -= microsecs();
		
	}
	void StopTimer( int timer )
	{
		timers[timer] += microsecs();
	}
	double GetTimer( int timer )
	{
		return timers[timer] * 1e-6;
	}
};
#endif
