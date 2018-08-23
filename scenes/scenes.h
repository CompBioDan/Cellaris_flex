// Simulation setup and data storing class

#ifndef SCENE_H_
#define SCENE_H_

#include <vector>
#include <climits>
#include "../utilities/datastore.h"

class Cell;

class Scene
{

protected: 

	double mDt; /** scene timestep */

	double mEndTime; /** scene end time */
	 
	unsigned mNumBirths; /** number of births during scene */

	int numActiveParticles; /** number of active particles in scene (used for FleX) */

	unsigned mSamplingTimestepMultiple; /** sampling (output) timestep ratio */

	std::string mSceneOutputDirectory; /** scene output directory */

public: 

	static Scene* Instance(); /** Scene instance */

	Scene(); /** Scene constructor */

	void addCell(Cell* cell); /** Method for adding a new cell to the scene */

	Cell* getCell(unsigned int cellID); /** Method to return specific cell entry */

	int getNumberCells(); /** return number of cells in scene */

	double getDt(); /** return the simulation timestep */

	unsigned getNumBirths(); /** return number of births in simulation */

	void setBirths(int births); /** set the number of births during simulation*/
	
	void setDt(double dt); /** set the timestep */

	void setEndTime(double endTime); /** set simulation end time*/

	void setOutputDirectory(std::string outputDirectory); /** set output directory */

	void setSamplingTimeStepMultiple(unsigned samplingTimestepMultiple); /** set the sampling (output) timestep ratio */

	void setNumberActiveParticles(int numactiveparticles); /** set number of active particles in the scene */

	int getNumberActiveParticles(); /** return the number of active particles in the scene */

	void Solve(); /** main solve method */

	static void Destroy(); /** destroy instance of the scene */

private:

	static Scene* mpInstance; /** pointer to the singleton instance of this class */

	CellarisVector<Cell*> cell_population; /** vector containing the population of cell pointers*/

};

#endif /*SCENE_H_*/