#include <cmath>
#include <iostream>
#include <fstream>

#include "scenes.h"
#include "../utilities/scenetime.h"
#include "../cells/cellpopulation.h"

// pointer to single instance
Scene* Scene::mpInstance = nullptr;


// @return a pointer to the scene object
Scene* Scene::Instance()
{
	if (mpInstance == nullptr)
	{
		mpInstance = new Scene;
		std::atexit(Destroy);
	}
	return mpInstance;
}

Scene::Scene() {}

double Scene::getDt()
{
	return mDt;
}

void Scene::setBirths(int births)
{
	mNumBirths = births;
}

unsigned Scene::getNumBirths()
{
	return mNumBirths;
}

void Scene::setDt(double dt)
{
	mDt = dt;
}

void Scene::setEndTime(double endTime)
{
	mEndTime = endTime;
}

void Scene::setOutputDirectory(std::string outputDirectory)
{
	mSceneOutputDirectory = outputDirectory;
}

void Scene::setSamplingTimeStepMultiple(unsigned samplingTimestepOutput)
{
	mSamplingTimestepMultiple = samplingTimestepOutput;
}

void Scene::Destroy()
{
	if (mpInstance)
	{
		delete mpInstance;
		mpInstance = nullptr;
	}
}

int Scene::getNumberCells()
{
	return cell_population.getCount();
}

int Scene::getNumberActiveParticles()
{
	return numActiveParticles;
}

void Scene::setNumberActiveParticles(int activeparticles)
{
	numActiveParticles = activeparticles;
}

void Scene::Solve()
{

	//// set up simulation time
	SceneTime* p_scene_time = SceneTime::Instance();
	double current_time = p_scene_time->GetTime();

	// set the end time and number of timesteps
	unsigned num_time_steps = (unsigned)((mEndTime - current_time) / mDt + 0.5);
	std::cout <<"tsteps "<< mDt << '\n';
	if (current_time > 0)
	{
		p_scene_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
	}
	else
	{
		p_scene_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
	}

	// set the current time
	double time_now = p_scene_time->GetTime();

	// correctly age all cells before main solve loop
	for (int i = 0; i < getNumberCells(); ++i)
	{
		cell_population[i]->readyToDivide();
	}
	std::cout << "end time: " << p_scene_time->GetEndTime() << '\n';
	std::cout << "current time: " << p_scene_time->GetTime() << '\n';
	std::cout << "at end? " << p_scene_time->HasFinished() << '\n';

	// Main solve loop
	while (!(p_scene_time->HasFinished()))
	{

		//UpdateCellPopulation();
		for (int i = 0; i < getNumberCells(); ++i)
		{

			if (cell_population[i]->readyToDivide())
			{
				addCell(cell_population[i]->divideCell());
			}
		}

		p_scene_time->IncrementTimeOneStep();

	}


}

void Scene::addCell(Cell* cell)
{
	cell_population.add(cell);
}

Cell* Scene::getCell(unsigned int cell_id)
{
	if (cell_id >= cell_population.getSize())
		return NULL; 
	else 
		return cell_population[cell_id];
}

