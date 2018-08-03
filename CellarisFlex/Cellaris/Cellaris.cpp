/*
	Initial Cellaris framework development code. 

	Cellaris is a framework incorporating cellular and subcellular modelling methods with the Nvidia FleX physics 
	library for particle-constaint based agent models

	Developments:
	1. Create our 'SceneTime' singleton for globally consistent time
		TESTS: A) create instance, b) set start time, c) set up time step, d) get time step, e) increment time step
	2. Create 'Cell' class as a container for all biological information about a cell
		TESTS: a) test cell constructor, b) test cell aging, c) test cell division, 
	


*/

#include <iostream>
#include "stdafx.h"
#include <cassert>
#define DEBUG

#include "cells/cell.h"
#include "utilities/scenetime.h"
#include "scenes/scenes.h"
//#include "utilities/timestepper.h"


//#include "cxxtest/cxxtest/TestSuite.h"
//#include "../cxxtest/TestSuite.h"

int main ()
{

	// Start of the testing for the Cellaris simulating framework
	std::cout << "Starting testing for the Cellaris framework..." << '\n' << '\n';

	////////////////// SceneTime testing 
	// Create scene time object - set sim length and number of time steps
	SceneTime* p_scene_time = SceneTime::Instance();

	// Set start time & test that it is allocated correctly
	double startTime = 0.0;
	std::cout << "Test 1: Setting scene start time." << '\n';
	p_scene_time->SetStartTime(startTime); // Set simulation start time to zero
	std::cout << "Set scene start time is: " << p_scene_time->GetStartTime() << " - Actual start time is: " << startTime << '\n' << '\n';
	
	// Set end time and timesteps & test allocation
	std::cout << "Test 2: Setting the scene end time and the timestep." << '\n';
	p_scene_time->SetEndTimeAndNumberOfTimeSteps(10.0, 3); // Set end time and number of timesteps
	std::cout << "Scene end time is: " << p_scene_time->GetEndTime() << " - Timestep is: " << p_scene_time->GetTimeStep() << '\n' << '\n';

	// Increment time step & test timestep incrementing
	std::cout << "Test 3: Checking time step incrementing." << '\n';
	std::cout << "Number of timesteps taken = " << p_scene_time->GetTimeStepsElapsed() << " - Current time is = " << p_scene_time->GetTime() << '\n';
	p_scene_time->IncrementTimeOneStep();
	std::cout << "Number of timesteps taken = " << p_scene_time->GetTimeStepsElapsed() << " - Current time is = " << p_scene_time->GetTime() << '\n' << '\n';

	// Destroy the scene time object - avoid memory leak
	SceneTime::Destroy();

	///////////////// Cell class testing
	// Create a cell
	Cell* newCell = new Cell();
	
	// Set cell data (birthtime, position and id)
	doubleVec3d newPos;
	newPos.x = 1.0; newPos.y = 1.0; newPos.z = 1.0;

	newCell->setBirthTime(-0.5); // set cell birth time
	newCell->setCellPos(newPos); // set cell position
	newCell->setCellId(0); // set cell id

	std::cout << "Test 4: Setting cell parameters" << '\n' << "Cell birth time = " << newCell->getBirthTime() <<  '\n';
	doubleVec3d retrieved_pos = newCell->getCellPos();
	std::cout << "Cell position = x: " << retrieved_pos.x << " y: " << retrieved_pos.y << " z: " << retrieved_pos.z << '\n';
	std::cout << "Cell ID = " << newCell->getCellId() << '\n' << '\n';

	// Test cell aging and timestep incrementing
	SceneTime* p_scene_time1 = SceneTime::Instance(); // Create new scenetime instance

	// Set start time & test
	p_scene_time1->SetStartTime(0.0); // Set simulation start time to zero

	// Set end time and timesteps & test
	double endTime = 10.0;
	unsigned int timeSteps = 20;
	p_scene_time1->SetEndTimeAndNumberOfTimeSteps(endTime, timeSteps); // Set end time and number of timesteps

	std::cout << "Test 5: Testing cell aging and timestep incrementing: " << '\n';
	std::cout << "Timestep size = " << endTime/timeSteps <<'\n';
	std::cout << "Cell age = " << newCell->getAge() << '\n';
	p_scene_time1->IncrementTimeOneStep();
	std::cout << "Increment 1 time step" << '\n' << "Cell age = " << newCell->getAge() << '\n';
	p_scene_time1->IncrementTimeOneStep();
	p_scene_time1->IncrementTimeOneStep();
	std::cout << "Increment 2 time steps" << '\n' << "Cell age = " << newCell->getAge() << '\n' << '\n';

	SceneTime::Destroy(); // Destroy scenetime instance to avoid memory leak

	// Test cell division
	Cell* newCell2 = new Cell();

	// Setup cell parameters (position, birthtime, and id)
	doubleVec3d nPos2;
	nPos2.x = 1.0; nPos2.y = 1.0; nPos2.z = 1.0;

	newCell2->setBirthTime(0.0); // set cell birth time
	newCell2->setCellPos(nPos2); // set cell position
	newCell2->setCellId(0); // set cell id

	SceneTime* p_scene_time2 = SceneTime::Instance(); // New scenetime instance
	p_scene_time2->SetStartTime(0.0); // Set simulation start time to zero
	p_scene_time2->SetEndTimeAndNumberOfTimeSteps(15.0, 3); // Set end time and number of timesteps

	Scene* p_scene = Scene::Instance(); // New scene instance
	p_scene->setBirths(0);

	std::cout << "Test 6: Test the readyToDivide() method (set end of cell-cycle to 10h, timestep is 5h)" << '\n';
	p_scene_time2->IncrementTimeOneStep();
	std::cout << "Ready to divide?: " << newCell2->readyToDivide() << " cell age = " << newCell2->getAge() << '\n';
	p_scene_time2->IncrementTimeOneStep();
	std::cout << "Ready to divide?: " << newCell2->readyToDivide() << " cell age = " << newCell2->getAge() << '\n' << '\n';

	std::cout << "Test 7: Test the scene tracking of cell birth events. " << '\n';
	std::cout << "Number of births: " << p_scene->getNumBirths() << '\n';
	Cell* p_daughter_cell = newCell2->divideCell();
	std::cout << "Cell division occured.." << '\n';
	std::cout << "Number of births: " << p_scene->getNumBirths() << '\n' << '\n';
	
	// Test scene ability to track ages and id
	std::cout << "Test 8: Test age and id tracking for the cells. " << '\n';
	std::cout << "Mother cell age: " << newCell2->getAge() << '\n';
	std::cout << "Mother cell id: " << newCell2->getCellId() << '\n';
	std::cout << "Daughter cell age: " << p_daughter_cell->getAge() << '\n';
	std::cout << "Daughter cell id: " << p_daughter_cell->getCellId() << '\n' << '\n';

	// Destroy scenetime and scene instances
	SceneTime::Destroy();
	Scene::Destroy();

	// Testing a multiple cell division example
	double end_time = 25.0;
	int time_steps = 25;

	// Create new scene time instance and set initial time, end time and timesteps
	SceneTime* p_scene_time3 = SceneTime::Instance();
	p_scene_time3->SetStartTime(0.0);
	p_scene_time3->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);

	// Create initial cell and set birthtime, cellid and a position
	Cell* p_cell(new Cell());
	p_cell->setBirthTime(p_scene_time3->GetTime());
	p_cell->setCellId(0);
	p_cell->setCellPos(nPos2);

	// Create a vector of cells, newly born cells and track times
	std::vector<Cell*> cells;
	std::vector<Cell*> newly_born;
	std::vector<double> times(time_steps);

	// Add intial cell to the cell vector
	cells.push_back(p_cell);
	std::vector<Cell*>::iterator cell_iterator;

	// Create our scene instance and set initial birth count to zero
	Scene* p_scene2 = Scene::Instance();
	p_scene2->setBirths(0);

	int i = 0;
	std::cout << "Test 9: Testing multiple division in a simulation loop." << '\n';
	while (p_scene_time3->GetTime() < end_time)
	{
		p_scene_time3->IncrementTimeOneStep(); // Increment simulation
		cell_iterator = cells.begin(); 

		// Loop through all cells and check whether cell is ready to divide, if so then create daughter cell
		while (cell_iterator < cells.end())
		{
			if ((*cell_iterator)->readyToDivide())
			{
				newly_born.push_back((*cell_iterator)->divideCell());
			}
			++cell_iterator;
		}

		// Add the newly born cells to the cells vector
		cell_iterator = newly_born.begin();
		while (cell_iterator < newly_born.end())
		{
			cells.push_back(*cell_iterator);
			++cell_iterator;
		}
		newly_born.clear();

		// Track simulation time
		times[i] = p_scene_time3->GetTime();
		i++;

		std::cout << "Time " << p_scene_time3->GetTime() << ": Cell count: " << cells.size() <<  '\n';

	}

	std::cout << '\n'<< "Comparing true births with scene tracking of births:" << '\n';
	std::cout << "True birth count = " << cells.size()-1 << '\n';
	std::cout << "Scene birth count = " << p_scene->getNumBirths() << '\n' << '\n';
	std::cout << "Id of mother cell = " << cells[0]->getCellId() << '\n';
	std::cout << "Id of daughter cell = " << cells[1]->getCellId() << '\n';
	std::cout << "Id of granddaughter cell = " << cells[2]->getCellId() << '\n';

	SceneTime::Destroy();
	Scene::Destroy();
	std::cin.ignore();
	return 0;
}





