
#include "cellpopulation.h"


CellPopulation::CellPopulation(std::vector<Cell*>& rCells)
	: mCells(rCells.begin(), rCells.end())
{
	std::vector<Cell*>().swap(rCells);
}

void CellPopulation::InitialiseCells()
{
	// Can set up the initializing of cell-cycle, cell positions etc here
}

std::list<Cell*>& CellPopulation::rGetCells()
{
	return mCells;
}

unsigned CellPopulation::getNumAllCells()
{
	return mCells.size();
}



