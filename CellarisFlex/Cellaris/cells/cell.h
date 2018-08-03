// Generic cell class to contain the representation of 'cells' constructed from particles and particle-particle constraints
//
// Base cell class contains information about the cells particle construction which is passed to FleX to solve physical interactions
//
// 
//

#ifndef CELL_H_
#define CELL_H_

#include "../stdafx.h"
//#include "../core/vec3.h"
#include "../utilities/datastore.h"
#include "../utilities/scenetime.h"
#include "../scenes//scenes.h"

#include <vector>
#include <cmath>

class Cell {

private:

	bool mEndofCycle; // has the cell reached the end of its cell-cycle -> check to pass to initiate cell division

public:

	Cell(); // create a new cell

	//virtual ~Cell();

	//virtual void update();

	void setBirthTime(double birthTime); // set the cells birthtime

	double getBirthTime() const;

	double getAge() const;

	void setCellPos(doubleVec3d position);

	doubleVec3d getCellPos() const;

	bool readyToDivide();

	void setCellId(unsigned cellId);

	void resetCell();

	unsigned getCellId() const;

	Cell* divideCell();

	// Cell data
	double mCellBirthTime;

	std::vector<int> mIndices; // list of indices for particles belonging to cell

	doubleVec3d mPosition; // position of first particle in cell (arbitrary)

	unsigned mCellId; // identifier id for the cell

	double mCellAge; // age of the cell


};

#endif /* CELL_H_ */