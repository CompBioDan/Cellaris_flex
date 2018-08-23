#pragma once
#pragma once


// TEST SCENE FOR SWARM CELL GROWTH

#include <random>
class SwarmBacteriaGrowth : public Scene
{

public:
	SwarmBacteriaGrowth(const char* name) : Scene(name) {}

	struct Bacteria
	{
		// particles
		int particle_offset; // index of first particle in cell
		int number_particles; // number of particles in the cell

							  // springs
		int spring_offset; // index of first spring length/stiffness value
		int number_springs; // number of springs within cell

							// cell characteristics
		int cell_type; // cell type flag
		double cell_age; // cell age
		int grow; // cell growth flag
		int divide; // cell divide flag
		int cell_id; // cell id value
	};

	virtual void Initialize()
	{

		// Scene domain setup parameters///////////////////////
		float min_X = 0.0f; float min_Y = 0.0f; float min_Z = 0.0f;
		//float max_X = 40.0f; float max_Y = 20.0f; float max_Z = 40.0f;
		float max_X = 20.0f; float max_Y = 10.0f; float max_Z = 20.0f;

		g_sceneLower = Vec3(min_X, min_Y, min_Z);
		g_sceneUpper = Vec3(max_X, max_Y, max_Z);

		g_camPos = Vec3(10.0f, 20.0f, 10.0f);
		g_camAngle = Vec3(-DegToRad(30.0f), -DegToRad(27.0f), 0.0f);
		///////////////////////////////////////////////////////

		// Limit on max number of particles/springs in the simulation
		g_numExtraParticles = 1000000;
		int g_numExtraSprings = 30;
		/////////////////////////////////////////////////////////////

		//  Bacteria parameters//////////////////////////////////////
		radius = 0.2f; // Radius of 'particles' forming the cells
		length = 0.4f; // Length of the 'rod' cells
		numbersegments = 5;
		growthProb = 0.01;
		/////////////////////////////////////////////////////////////

		// Simulation parameters/////////////////////////////////////
		g_numSubsteps = 3;
		g_params.numIterations = 5;
		g_params.radius = radius;
		g_params.dynamicFriction = 0.2f;
		g_params.dissipation = 0.001f;
		g_params.sleepThreshold = g_params.radius*0.2f;
		g_params.relaxationFactor = 1.0f;
		g_params.restitution = 0.0f;
		g_params.shapeCollisionMargin = 0.01f;
		g_params.particleCollisionMargin = g_params.radius*0.05f;
		g_params.gravity[1] = 0.0f;
		g_params.numPlanes = 6;
		g_lightDistance *= 0.5f;
		g_drawPoints = true;
		g_drawRopes = true;
		g_drawSprings = true;
		dt = 1 / (double)60;
		mFrame = 0;


		// Cell population parameters
		cellCount = 1000; // starting cell count
		numSprings = 0; // number of springs in the scene
		maxCellParticles = 6; // Max number of particles in each cell
		numExtraParticles = 0; // Number of extra particles to add to the simulation
		sceneParticleCount = 0; // Number of particles in the scene (used for initialising buffer tranfer holder
		/////////////////////////////////////////////////////////////

		// Create the initial cells//////////////////////////////////
		int offsetCounter = 0; // Loop variable for setting first buffer index of cell particle data
		int number_springser = 0; // Loop variable for setting the first buffer index of cell spring data
		numbersegments = 2;

		// Create the cells 
		for (int i = 0; i < cellCount; ++i)
		{
			Bacteria bac;
			bac.particle_offset = offsetCounter;
			bac.spring_offset = number_springser;
			bac.number_particles = int(ceil(length / radius)) + 1;
			bac.number_springs = 0;
			bac.cell_type = 1;
			bac.cell_age = 0.0;
			bac.cell_id = i;
			bac.grow = 0;
			bac.divide = 0;

			CreateBacteria(bac, Vec3(Randf(min_X + length, max_X - length), Randf(min_Y + length, max_Y - length), Randf(min_Z + length, max_Z - length)), Vec3(1.0f, 0.0f, 0.0f), 0.25f, int(ceil(length / radius)), length, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide));
			g_bac.push_back(bac);

			offsetCounter += bac.number_particles;
			number_springser += bac.number_springs;

			sceneParticleCount += bac.number_particles;
		}
	}

	virtual void Sync()
	{
		// update solver
		NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
	}


	void Update()
	{
		// Age the cells in the scene
		for (int cell_counter = 0; cell_counter < cellCount; cell_counter++)
		{
			g_bac[cell_counter].cell_age += dt;
		}

		// Flag cells for growth
		int flagged_grow = 0;
		for (int cell_counter = 0; cell_counter < cellCount; cell_counter++)
		{
			float grow_prob = Randf(0.0f, 1.0f);
			if (grow_prob < growthProb & g_bac[cell_counter].number_particles < maxCellParticles)
			{
				g_bac[cell_counter].grow = 1;
				flagged_grow += 1;
			}
		}

		// Grow the cells flagged for growth
		for (int i = 0; i < cellCount; i++)
		{
			if (g_bac[i].grow == 1)
			{
				GrowCells(i);
			}
		}

		++mFrame;
	}

	void GrowCells(int cellID)
	{
		//int predictedParticleCount = 

		int cellFirstParticleIndex;
		int cellFirstSpringIndex;

		// Resize the interim buffers to zero to be rebuilt
		numParticleGrowth = 0;
		particles.resize(0);
		velocities.resize(0);
		phases.resize(0);

		springstiffness.resize(0);
		springlengths.resize(0);
		springindices.resize(0);

		// Looping variable for rebuilding the active particle and spring buffers
		numActiveParticles = 0;
		numActiveSprings = 0;

		// Copy cells into relevant holders from the buffers
		for (int cell_counter = 0; cell_counter < cellCount; cell_counter++)
		{
			// Index of first particle in cell
			cellFirstParticleIndex = g_bac[cell_counter].particle_offset - numParticleGrowth;
			cellFirstSpringIndex = g_bac[cell_counter].spring_offset - (2 * numParticleGrowth);

			int numberParticlesInCell = g_bac[cell_counter].number_particles;
			int numberSpringsInCell = g_bac[cell_counter].number_springs;
			int endParticleIndex = (cellFirstParticleIndex + numberParticlesInCell) - 1;

			// Create temporary buffers for the particles in the simulation
			for (int j = 0; j < numberParticlesInCell; ++j)
			{
				particles.push_back(g_buffers->positions[cellFirstParticleIndex + j]);
				velocities.push_back(g_buffers->velocities[cellFirstParticleIndex + j]);
				phases.push_back(g_buffers->phases[cellFirstParticleIndex + j]);
			}

			for (int j = 0; j < numberSpringsInCell; ++j)
			{
				springstiffness.push_back(g_buffers->springStiffness[cellFirstSpringIndex + j]);
				springlengths.push_back(g_buffers->springLengths[cellFirstSpringIndex + j]);
			}

			for (int j = 0; j < numberSpringsInCell * 2; ++j)
			{
				springindices.push_back(g_buffers->springIndices[cellFirstSpringIndex * 2 + j] + numParticleGrowth);
			}


			numActiveParticles += numberParticlesInCell;
			numActiveSprings += numberSpringsInCell;

			particleSpacing.x = g_buffers->positions[endParticleIndex].x - g_buffers->positions[endParticleIndex - 1].x;
			particleSpacing.y = g_buffers->positions[endParticleIndex].y - g_buffers->positions[endParticleIndex - 1].y;
			particleSpacing.z = g_buffers->positions[endParticleIndex].z - g_buffers->positions[endParticleIndex - 1].z;

			// If cell is flagged for growth then add a new particle to the end of the bacteria and add to the buffer holders
			if (cell_counter == cellID)
			{
				// Set the particle to end of the current cell 
				particles.push_back(Vec4(g_buffers->positions[endParticleIndex].x + (particleSpacing.x), g_buffers->positions[endParticleIndex].y + (particleSpacing.y), g_buffers->positions[endParticleIndex].z + (particleSpacing.z), 1.0f));
				velocities.push_back(g_buffers->velocities[endParticleIndex]);
				phases.push_back(g_buffers->phases[endParticleIndex]);

				// spring 1 (new to end)
				springstiffness.push_back(0.25f);
				springlengths.push_back(Length(Vec3(g_buffers->positions[endParticleIndex]) - Vec3(particles[numActiveParticles])));
				springindices.push_back(numParticleGrowth + (endParticleIndex));
				springindices.push_back(numActiveParticles);

				// spring 2 (new to end -1)
				springstiffness.push_back(0.25f);
				springlengths.push_back(Length(Vec3(g_buffers->positions[endParticleIndex - 1]) - Vec3(particles[numActiveParticles])));
				springindices.push_back(numParticleGrowth + (endParticleIndex - 1));
				springindices.push_back(numActiveParticles);

				// Increment cell particle and spring counters
				g_bac[cell_counter].number_particles += 1;
				g_bac[cell_counter].number_springs += 2;

				for (int cc = cell_counter + 1; cc < cellCount; cc++)
				{
					g_bac[cc].particle_offset += 1;
					g_bac[cc].spring_offset += 2;
				}

				numActiveParticles += 1;
				numActiveSprings += 2;
				grow = 0;
				numParticleGrowth += 1;

			}

			g_bac[cell_counter].grow = 0;

		}

		// Assign the holders to the FleX buffers
		g_buffers->positions.assign(&particles[0], particles.size());
		g_buffers->velocities.assign(&velocities[0], velocities.size());
		g_buffers->phases.assign(&phases[0], phases.size());

		g_buffers->springIndices.assign(&springindices[0], springindices.size());
		g_buffers->springLengths.assign(&springlengths[0], springlengths.size());
		g_buffers->springStiffness.assign(&springstiffness[0], springstiffness.size());

		// Extend the active indices to include any new added particles
		g_buffers->activeIndices.resize(numActiveParticles);
		for (int i = 0; i < numActiveParticles; ++i)
			g_buffers->activeIndices[i] = i;

	}

	void ForceRestart()
	{

		int offsetCounter = 0; // Loop variable for setting first buffer index of cell particle data
		int number_springser = 0; // Loop variable for setting the first buffer index of cell spring data
		sceneParticleCount = 0;
		for (int i = 0; i < cellCount; ++i)
		{
			g_bac[i].particle_offset = offsetCounter;
			g_bac[i].spring_offset = number_springser;
			g_bac[i].number_particles = int(ceil(length / radius)) + 1;
			g_bac[i].cell_type = 1;
			g_bac[i].cell_age = 0.0;
			g_bac[i].number_springs = (int(ceil(length / radius)) * 2) - 3;
			g_bac[i].cell_id = i;
			g_bac[i].grow = 0;
			g_bac[i].divide = 0;

			offsetCounter += g_bac[i].number_particles;
			number_springser += g_bac[i].number_springs;

			sceneParticleCount += g_bac[i].number_particles;

		}

		//g_bac[1].divide = 1;

	}

	void CreateBacteria(Bacteria& bac, Vec3 start, Vec3 dir, float stiffness, int segments, float length, int phase, float spiralAngle = 0.0f, float invmass = 1.0f, float give = 0.075f)
	{

		//bac.mIndices.push_back(int(g_buffers->positions.size()));

		g_buffers->positions.push_back(Vec4(start.x, start.y, start.z, invmass));
		g_buffers->velocities.push_back(0.0f);
		g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));

		Vec3 left, right;
		BasisFromVector(dir, &left, &right);

		float segmentLength = length / segments;
		Vec3 spiralAxis = dir;
		float spiralHeight = spiralAngle / (2.0f*kPi)*(length / segments);

		if (spiralAngle > 0.0f)
			dir = left;

		Vec3 p = start;

		for (int i = 0; i < segments; ++i)
		{
			int prev = int(g_buffers->positions.size()) - 1;

			p += dir * segmentLength;

			// rotate 
			if (spiralAngle > 0.0f)
			{
				p += spiralAxis * spiralHeight;

				dir = RotationMatrix(spiralAngle, spiralAxis)*dir;
			}

			//bac.mIndices.push_back(int(g_buffers->positions.size()));

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);

			// stretch spring
			CreateSpring(prev, prev + 1, stiffness, give);
			bac.number_springs += 1;
			numSprings += 1;

			// bending spring
			if (i > 0)
			{
				CreateSpring(prev - 1, prev + 1, stiffness*0.5f, give);
				bac.number_springs += 1;
				numSprings += 1;
			}
		}

		bac.number_particles = segments + 1;

	}


	vector<Bacteria> g_bac; // vector containing all of the bacteria cells

	Vec3 particleSpacing;
	int numActiveParticles; // integer counter for rebuilding the buffers in the grow step
	int numActiveSprings; // integer counter for rebuilding the buffers in the grow step

	int mFrame; // tracking frame number
	double dt;
	double time;
	int phase; // tracking cell phase for creating new particles 
	float velocity; // tracking velocity for creating new particles
	int sceneParticleCount;
	int numSprings; // number of springs in sim
	int cellCount; // number of cells in sim

	int grow; // growth flag

	int maxCellParticles; // maximum number of particles in cell

	int numExtraParticles;

	// Interim buffers for the cells
	std::vector<Vec4> particles;
	std::vector<Vec3> velocities;
	std::vector<int> phases;
	std::vector<float> springstiffness;
	std::vector<float> springlengths;
	std::vector<int> springindices;

	int numParticleGrowth; // number of particles added in Update() function
	int numbersegments;

	float length;
	float radius;
	double growthProb;
};