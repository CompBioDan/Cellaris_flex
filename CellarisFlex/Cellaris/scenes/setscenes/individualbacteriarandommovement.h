#include <random>
class IndividualBacteriaRandomMovement : public Scene
{

public:
	IndividualBacteriaRandomMovement(const char* name) : Scene(name) {}

	struct Bacteria
	{
		// particles
		int particle_offset; // index of first particle in cell
		int number_particles; // number of particles in the cell

							  // springs
		int spring_offset;
		int number_springs; // number of springs within cell

							// cell characteristics
		int cell_type; // cell type flag
		double cell_age; // cell age
		int grow;
		int divide;
		int cell_id;

		float cell_velocity;
	};

	virtual void Initialize()
	{

		// Scene domain setup parameters///////////////////////
		float min_X = 0.0f; float min_Y = 0.0f; float min_Z = 0.0f;
		//float max_X = 40.0f; float max_Y = 20.0f; float max_Z = 40.0f;
		float max_X = 10.0f; float max_Y = 5.0f; float max_Z = 10.0f;

		g_sceneLower = Vec3(min_X, min_Y, min_Z);
		g_sceneUpper = Vec3(max_X, max_Y, max_Z);
		///////////////////////////////////////////////////////

		// Limit on max number of particles/springs in the simulation
		g_numExtraParticles = 1000000;
		int g_numExtraSprings = 30;
		/////////////////////////////////////////////////////////////

		//  Bacteria parameters//////////////////////////////////////
		float radius = 0.2f; // Radius of 'particles' forming the cells
		float length = 0.4f; // Length of the 'rod' cells
		numbersegments = 5;
		mVelocity = 5.0f;
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
		cellCount = 1; // starting cell count
		numSprings = 0; // number of springs in the scene
		maxCellParticles = 6; // Max number of particles in each cell
		numExtraParticles = 0; // Number of extra particles to add to the simulation
		sceneParticleCount = 0; // Number of particles in the scene (used for initialising buffer tranfer holder
								/////////////////////////////////////////////////////////////

		// Create the initial cells//////////////////////////////////
		int offsetCounter = 0; // Loop variable for setting first buffer index of cell particle data
		int number_springser = 0; // Loop variable for setting the first buffer index of cell spring data
		numbersegments = 2;
		for (int i = 0; i < cellCount; ++i)
		{
			Bacteria bac;
			bac.particle_offset = offsetCounter;
			bac.spring_offset = number_springser;
			bac.number_particles = int(ceil(length / radius))+1;
			bac.number_springs = 0;
			bac.cell_type = 1;
			bac.cell_age = 0.0;
			bac.cell_id = i;
			bac.grow = 0;
			bac.divide = 0;
			bac.cell_velocity = 5.0f;

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

	virtual void DoGui() // Include a parameter varying GUI
	{
		if (imguiSlider("Bacteria Velocity", &mVelocity, 0.0f, 20.0f, 0.1f))
		{
			for (int i = 0; i < int(g_bac.size()); ++i)
				g_bac[i].cell_velocity = mVelocity;
		}
	}

	void Update()
	{

		// Force cell-reset
		if (mFrame == 0) {
			ForceRestart();
		}

		// Age the cells in the scene
		for (int cell_counter = 0; cell_counter < cellCount; cell_counter++) {
			g_bac[cell_counter].cell_age += dt;
		}

		for (int i = 0; i < cellCount; i++) {
			MoveCells(i);
		}

		++mFrame;
	}

	void MoveCells(int cell_id)
	{
		//// 1. Swim
		// 1.1 Find vector in forward direction
		// set number of particles in cell
		int cellFirstParticleIndex;
		int cellFirstSpringIndex;

		// Index of first particle in cell
		cellFirstParticleIndex = g_bac[cell_id].particle_offset;
		cellFirstSpringIndex = g_bac[cell_id].spring_offset;

		int numberParticlesInCell = g_bac[cell_id].number_particles;
		int numberSpringsInCell = g_bac[cell_id].number_springs;
		int endParticleIndex = (cellFirstParticleIndex + numberParticlesInCell) - 1;

		//int particleCount = g_bac[cell_id].number_particles;
		//int particle_offset = g_bac[cell_id].particle_offset + particleCount;
		Vec2 tumblePoints = Vec2(cellFirstParticleIndex, endParticleIndex);

		// calculate the average direction vector of the cell
		Vec3 dirvector1 = Vec3(g_buffers->positions[endParticleIndex]) - Vec3(g_buffers->positions[endParticleIndex - 1]);
		Vec3 dirvector2 = Vec3(g_buffers->positions[endParticleIndex - 1]) - Vec3(g_buffers->positions[endParticleIndex - 2]);
		Vec3 averageDirectionVector;
		averageDirectionVector.x = (dirvector1.x + dirvector2.x) / 2;
		averageDirectionVector.y = (dirvector1.y + dirvector2.y) / 2;
		averageDirectionVector.z = (dirvector1.z + dirvector2.z) / 2;
		// Pick a velocity to apply from a uniform distribution from 0 to a maxVel)
		float ins_vel = Randf(0.0f, g_bac[cell_id].cell_velocity);
		float tumble_vel = Randf(0.0f, 20.0f);

		//// 2. Tumble
		//
		int randomTumbleDirection = floor(Randf(0, 6));
		float tumbleProbability = Randf();

		if (tumbleProbability > 0.05) {
			for (int i = 0; i < numberParticlesInCell; ++i) {
				g_buffers->velocities[g_bac[cell_id].particle_offset + i].x = ins_vel * averageDirectionVector.x;
				g_buffers->velocities[g_bac[cell_id].particle_offset + i].y = ins_vel * averageDirectionVector.y;
				g_buffers->velocities[g_bac[cell_id].particle_offset + i].z = ins_vel * averageDirectionVector.z;
			}
		}
		else {
			if (randomTumbleDirection == 0) {
				g_buffers->velocities[endParticleIndex].x = tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[endParticleIndex].y = 0;
				g_buffers->velocities[endParticleIndex].z = -tumble_vel * averageDirectionVector.x;

				g_buffers->velocities[tumblePoints[0]].x = -tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[tumblePoints[0]].y = 0;
				g_buffers->velocities[tumblePoints[0]].z = tumble_vel * averageDirectionVector.x;
			}
			else if (randomTumbleDirection == 1) {
				g_buffers->velocities[endParticleIndex].x = 0;
				g_buffers->velocities[endParticleIndex].y = tumble_vel * averageDirectionVector.z;
				g_buffers->velocities[endParticleIndex].z = -tumble_vel * averageDirectionVector.y;

				g_buffers->velocities[tumblePoints[0]].x = 0;
				g_buffers->velocities[tumblePoints[0]].y = -tumble_vel * averageDirectionVector.z;
				g_buffers->velocities[tumblePoints[0]].z = tumble_vel * averageDirectionVector.y;
			}
			else if (randomTumbleDirection == 2) {
				g_buffers->velocities[endParticleIndex].x = tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[endParticleIndex].y = -tumble_vel * averageDirectionVector.x;
				g_buffers->velocities[endParticleIndex].z = 0;

				g_buffers->velocities[tumblePoints[0]].x = -tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[tumblePoints[0]].y = tumble_vel * averageDirectionVector.x;
				g_buffers->velocities[tumblePoints[0]].z = 0;

			}
			else if (randomTumbleDirection == 3) {
				g_buffers->velocities[endParticleIndex].x = -tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[endParticleIndex].y = 0;
				g_buffers->velocities[endParticleIndex].z = tumble_vel * averageDirectionVector.x;

				g_buffers->velocities[tumblePoints[0]].x = tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[tumblePoints[0]].y = 0;
				g_buffers->velocities[tumblePoints[0]].z = -tumble_vel * averageDirectionVector.x;
			}
			else if (randomTumbleDirection == 4) {
				g_buffers->velocities[endParticleIndex].x = 0;
				g_buffers->velocities[endParticleIndex].y = -tumble_vel * averageDirectionVector.z;
				g_buffers->velocities[endParticleIndex].z = tumble_vel * averageDirectionVector.y;

				g_buffers->velocities[tumblePoints[0]].x = 0;
				g_buffers->velocities[tumblePoints[0]].y = tumble_vel * averageDirectionVector.z;
				g_buffers->velocities[tumblePoints[0]].z = -tumble_vel * averageDirectionVector.y;
			}
			else if (randomTumbleDirection == 5) {
				g_buffers->velocities[endParticleIndex].x = -tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[endParticleIndex].y = tumble_vel * averageDirectionVector.x;
				g_buffers->velocities[endParticleIndex].z = 0;

				g_buffers->velocities[tumblePoints[0]].x = tumble_vel * averageDirectionVector.y;
				g_buffers->velocities[tumblePoints[0]].y = -tumble_vel * averageDirectionVector.x;
				g_buffers->velocities[tumblePoints[0]].z = 0;

			}
		}

	}

	void ForceRestart()
	{

		//int offsetCounter = 0; // Loop variable for setting first buffer index of cell particle data
		//int number_springser = 0; // Loop variable for setting the first buffer index of cell spring data
		//sceneParticleCount = 0;
		//for (int i = 0; i < cellCount; ++i)
		//{

		//	g_bac[i].particle_offset = offsetCounter;
		//	g_bac[i].spring_offset = number_springser;
		//	g_bac[i].number_particles = numbersegments + 1;
		//	g_bac[i].cell_type = 1;
		//	g_bac[i].cell_age = 0.0;
		//	g_bac[i].number_springs = ((numbersegments + 1)*2)-3;
		//	g_bac[i].cell_id = i;
		//	g_bac[i].grow = 0;
		//	g_bac[i].divide = 1;

		//	offsetCounter += g_bac[i].number_particles;
		//	number_springser += g_bac[i].number_springs;

		//	sceneParticleCount += g_bac[i].number_particles;

		//}

		//g_bac[1].divide = 1;

	}

	void CreateBacteria(Bacteria& bac, Vec3 start, Vec3 dir, float stiffness, int segments, float length, int phase, float spiralAngle = 0.0f, float invmass = 1.0f, float give = 0.075f)
	{
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

			// stretch
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

	float mVelocity;
};