#include <cmath>
#include <iostream>
#include <fstream>

#include "scenes.h"

// pointer to single instance
Scene* Scene::mpInstance = nullptr;


// @return a pointer to the scene time object
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
