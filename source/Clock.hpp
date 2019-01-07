#ifndef CLOCK_HPP_
#define CLOCK_HPP_

#include <cinttypes>
#include <sys/time.h>

struct TimeInterval
{
	timespec begin;
	timespec end;
};

inline void beginTimeMeasurement(TimeInterval &ti)

{
	clock_gettime(CLOCK_MONOTONIC,&ti.begin);
}

inline void endTimeMeasurement(TimeInterval &ti)
{
	clock_gettime(CLOCK_MONOTONIC,&ti.end);
}

inline double timeIntervalToSeconds(const TimeInterval &ti)
{
	return (double)(ti.end.tv_sec-ti.begin.tv_sec)+static_cast<double>(ti.end.tv_nsec-ti.begin.tv_nsec)/1000000000.0;
}

#endif
