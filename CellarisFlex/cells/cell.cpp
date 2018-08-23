#include "../stdafx.h"
#include "cell.h"
#include <iostream>
#include <stdlib.h>
#include <random>

Cell::Cell() {


}

void Cell::setBirthTime(double birthTime) 
{
	cell_birth_time = birthTime;
}

double Cell::getAge() const 
{
	return SceneTime::Instance()->GetTime() - cell_birth_time;
}

double Cell::getBirthTime() const 
{
	return cell_birth_time;
}

//virtual void update();

doubleVec3d Cell::getCellPos() const
{
	return cell_position;
}

void Cell::setCellPos(doubleVec3d position) 
{
	cell_position = position;
}

void Cell::setCellId(unsigned cellId)
{
	cell_id = cellId;
}

unsigned Cell::getCellId() const
{
	return cell_id;
}

bool Cell::readyToDivide()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-1, 1);

	double cell_cycle_length = 10.0; /** TO CHANGE: new cell cycle class will allow for a cleaner allocation of a cell cycle length*/

	/** If cell has reached the end of its cell cycle then return true */
	if (Cell::getAge() >= (cell_cycle_length + dis(gen))) {
		return true;
	}
	else {
		return false;
	}
}

Cell* Cell::divideCell()
{
	/** Find number of births, used to allocate new ID for daughter cell */
	unsigned numberBirths = Scene::Instance()->getNumBirths();

	/** Reset the birthtime of the mother cell */
	resetCell();

	/** Create a new daughter cell */ 
	Cell* p_new_cell(new Cell());

	p_new_cell->setBirthTime(SceneTime::Instance()->GetTime()); /** Allocate the birth time of new cell to current simulation time*/
	p_new_cell->setCellId(numberBirths + 1); /** New daughter cell ID is set as the number of births */
	p_new_cell->setNumberParticles(3); /** TO CHANGE: number of particles will vary, depending on cell-type */
	p_new_cell->setParticleOffset(Scene::Instance()->getNumberActiveParticles()); /** Particle offset for flex buffers need seperate methods*/
	Scene::Instance()->setNumberActiveParticles(Scene::Instance()->getNumberActiveParticles() + 3); /** TO CHANGE: Number of active particles probably returned from Flex instance*/

	/** Increment the number of births during the scene */
	Scene::Instance()->setBirths(numberBirths + 1);

	/** Return the new cell */
	return p_new_cell;
}

void Cell::resetCell()
{
	/** reset the mother cell age following cell division */
	setBirthTime(SceneTime::Instance()->GetTime());
}

void Cell::setNumberParticles(int numparticles)
{
	number_particles = numparticles;
}

void Cell::setParticleOffset(int particleoffset)
{
	particle_offset = particleoffset;
}

int Cell::getParticleOffset()
{
	return particle_offset;
}