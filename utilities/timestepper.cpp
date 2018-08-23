#include <cassert>
#include <cmath>
#include <iostream>

#include "timestepper.h"

TimeStepper::TimeStepper(double startTime, double endTime, double dt)
	: mStart(startTime),
	mEnd(endTime),
	mDt(dt),
	mTotalTimeStepsTaken(0),
	mTime(startTime)
{

	try {
		// Output an error if start time is later than end time
		if (startTime > endTime)
		{
			throw 1;
		}
	}

	catch (int e)
	{
		std::cout << "ERROR: Simulation duration must be positive " << '\n';
	}


	mNextTime = CalculateNextTime();

}

double TimeStepper::CalculateNextTime()
{

	double next_time = mStart + (mTotalTimeStepsTaken + 1)*mDt;

	return next_time;

}

void TimeStepper::AdvanceOneTimeStep()
{

	mTotalTimeStepsTaken++;

	mTime = mNextTime;

	mNextTime = CalculateNextTime();

}

double TimeStepper::GetTime()
{
	return mTime;
}

double TimeStepper::GetNextTime()
{
	return mNextTime;
}

double TimeStepper::GetEndTime()
{
	return mEnd;
}

bool TimeStepper::IsTimeAtEnd() const
{
	return (mTime >= mEnd);
}

unsigned TimeStepper::GetTotalTimeStepsTaken() const
{
	return mTotalTimeStepsTaken;
}

double TimeStepper::GetNextTimeStep()
{
	double dt = mDt;

	if (mNextTime >= mEnd)
	{
		dt = mEnd - mTime;
	}

	return dt;
}

