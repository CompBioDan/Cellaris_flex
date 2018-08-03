#include "../stdafx.h"
#include "cell.h"
#include <iostream>

Cell::Cell() {}


void Cell::setBirthTime(double birthTime) 
{
	mCellBirthTime = birthTime;
}

double Cell::getAge() const 
{
	return SceneTime::Instance()->GetTime() - mCellBirthTime;
}

double Cell::getBirthTime() const 
{
	return mCellBirthTime;
}

//virtual void update();

doubleVec3d Cell::getCellPos() const
{
	return mPosition;
}

void Cell::setCellPos(doubleVec3d position) 
{
	mPosition = position;
}

void Cell::setCellId(unsigned cellId)
{
	mCellId = cellId;
}

unsigned Cell::getCellId() const
{
	return mCellId;
}

bool Cell::readyToDivide()
{
	double mCellAge = Cell::getAge();
	if (mCellAge >= 10){
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

	p_new_cell->setBirthTime(SceneTime::Instance()->GetTime());
	p_new_cell->setCellId(numberBirths + 1); /** New daughter cell ID is set as the number of births */

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
