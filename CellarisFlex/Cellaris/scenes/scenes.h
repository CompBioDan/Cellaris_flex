// Simulation setup and data storing class

#ifndef SCENE_H_
#define SCENE_H_

#include <vector>
#include <climits>

class Scene
{

protected: 

	// scene timestep
	double mDt; 

	// scene end time
	double mEndTime;

	// number of births during scene
	unsigned mNumBirths;

	// sampling (output) ratio 
	unsigned mSamplingTimestepMultiple;

	// simulation output directory
	std::string mSceneOutputDirectory;

public: 

	static Scene* Instance();

	// Constructor
	Scene(); // population object as input & initialise cellsdoublt 

	// return the simulation timestep
	double getDt();

	// return number of births in simulation
	unsigned getNumBirths();

	// set (for updating) number of births
	void setBirths(int births);

	// set the dt
	void setDt(double dt);

	// set sim end time
	void setEndTime(double endTime);

	// set output directory
	void setOutputDirectory(std::string outputDirectory);

	// set the sampling (output) timestep ratio 
	void setSamplingTimeStepMultiple(unsigned samplingTimestepMultiple);

	// Main solve method
	void Solve();

	static void Destroy();

private:

	// pointer to the singleton instance of this class
	static Scene* mpInstance;


};

#endif /*SCENE_H_*/