#include <cmath>
#include <cassert>
#include "scenetime.h"


// pointer to single instance
SceneTime* SceneTime::mpInstance = nullptr;

// shared pointed to the delegated class
std::shared_ptr<TimeStepper> SceneTime::mpTimeStepper;

// @return a pointer to the scene time object
SceneTime* SceneTime::Instance() 
{
	if (mpInstance == nullptr)
	{
		mpInstance = new SceneTime;
		mpTimeStepper.reset();
		std::atexit(Destroy);
	}
	return mpInstance;
}

// default scene time constructor
SceneTime::SceneTime()
	:
	mStartTime(0.0)
{
	// make sure there is only a single instance 
	assert(mpInstance == nullptr);
}

// destroys the current SceneTime instance - called before program exit to avoid a memory leak
void SceneTime::Destroy() 
{
	if (mpInstance)
	{
		delete mpInstance;
		mpInstance = nullptr;
	}
}

// get the simulation time step
double SceneTime::GetTimeStep() const
{
	return mpTimeStepper->GetNextTimeStep();
}

// increment simulation by a single time step
void SceneTime::IncrementTimeOneStep()
{
	mpTimeStepper->AdvanceOneTimeStep();
}

// get the number of time steps that have elapsed since start of the simulation
unsigned SceneTime::GetTimeStepsElapsed() const
{
	return mpTimeStepper->GetTotalTimeStepsTaken();
}

// get the simulation time
double SceneTime::GetTime() const
{
	if (mpTimeStepper)
	{
		return mpTimeStepper->GetTime();
	}
	return mStartTime;
}

// set the start time of the simulation
void SceneTime::SetStartTime(double startTime)
{
	mStartTime = startTime;
}

double SceneTime::GetStartTime() const
{
	return mStartTime;
}

double SceneTime::GetEndTime() const
{
	return mpTimeStepper->GetEndTime();
}

// set the end time of the simulation and number of timesteps
void SceneTime::SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation)
{
	mpTimeStepper.reset(new TimeStepper(mStartTime, endTime, (endTime - mStartTime) / totalTimeStepsInSimulation));
}

void SceneTime::ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun)
{
	/** start time is the timestepper current time*/
	mStartTime = mpTimeStepper->GetTime();

	assert(rEndTime > mStartTime);

	/** reset the machinery that works out the time */
	mpTimeStepper.reset(new TimeStepper(mStartTime, rEndTime, (rEndTime - mStartTime) / rNumberOfTimeStepsInThisRun));
}

bool SceneTime::HasFinished() const
{
	return(mpTimeStepper->IsTimeAtEnd());
}

