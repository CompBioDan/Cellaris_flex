/** 
*
* We define a 'Scene time' object to store the simulation time. 
* This object uses a singleton pattern to make sure we have a globally consistent time value.
*
* Start time, end time and number of time steps must be set before the time can be incremented and returned.
* 
* This class is used in conjunction with the IncrementTimeOneStep() and GetTime() calls.
*
 **/
#ifndef SCENETIME_H_
#define SCENETIME_H_

#include "timestepper.h"
#include <memory>

class SceneTime {

public: 

	/** 
	* @return a pointer to the scene time object
	*/
	static SceneTime* Instance();

	/** Sets the end time and number of timesteps - called after SetStartTime() but before other methods */
	void SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalSimulationTimeSteps);

	void ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun);

	/** @return the timestep for the simulation, set in earlier calls */
	double GetTimeStep() const;

	/** Increments the simulation time by a single time step */
	void IncrementTimeOneStep();

	/** @returns the number of timesteps that have elapsed in simulation thus far */
	unsigned GetTimeStepsElapsed() const;

	/** @returns the simulation time */
	double GetTime() const;

	/** Destroys the current SceneTime instance */
	static void Destroy();

	/** Set the simulation start time */
	void SetStartTime(double startTime);

	/** @return the simulation start time */
	double GetStartTime() const;

	/** @return the simulation end time */
	double GetEndTime() const;

	/** Has the simulation finished running */
	bool HasFinished() const;

protected:

	/** default scene time constructor -> must set start time, end time and number of timesteps before using object */
	SceneTime();

private:

	/** Pointer to the singleton instance of this class */
	static SceneTime* mpInstance;

	/** Delegate all time stepping to a TimeStepper class */
	static std::shared_ptr<TimeStepper> mpTimeStepper;

	/** Stores time at which simulation started */
	double mStartTime;


};

#endif /*SCENETIME_H_*/