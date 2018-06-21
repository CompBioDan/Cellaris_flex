#pragma once


class CellInflatable : public Scene
{
public:

	CellInflatable(const char* name) : Scene(name) {}

	virtual ~CellInflatable()
	{
		for (size_t i = 0; i < mCloths.size(); ++i)
			delete mCloths[i];
	}

	void AddCell(const Mesh* mesh, float overPressure, int phase)
	{
		const int startVertex = g_buffers->positions.size();

		// add mesh to system
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(mesh->m_positions[i]);

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);
		}

		int triOffset = g_buffers->triangles.size();
		int triCount = mesh->GetNumFaces();

		g_buffers->inflatableTriOffsets.push_back(triOffset / 3);
		g_buffers->inflatableTriCounts.push_back(mesh->GetNumFaces());
		g_buffers->inflatablePressures.push_back(overPressure);

		for (size_t i = 0; i < mesh->m_indices.size(); i += 3)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];

			Vec3 n = -Normalize(Cross(mesh->m_positions[b] - mesh->m_positions[a], mesh->m_positions[c] - mesh->m_positions[a]));
			g_buffers->triangleNormals.push_back(n);

			g_buffers->triangles.push_back(a + startVertex);
			g_buffers->triangles.push_back(b + startVertex);
			g_buffers->triangles.push_back(c + startVertex);
		}

		// create a cloth mesh using the global positions / indices
		ClothMesh* cloth = new ClothMesh(&g_buffers->positions[0], g_buffers->positions.size(), &g_buffers->triangles[triOffset], triCount * 3, 0.8f, 1.0f);

		for (size_t i = 0; i < cloth->mConstraintIndices.size(); ++i)
			g_buffers->springIndices.push_back(cloth->mConstraintIndices[i]);

		for (size_t i = 0; i < cloth->mConstraintCoefficients.size(); ++i)
			g_buffers->springStiffness.push_back(cloth->mConstraintCoefficients[i]);

		for (size_t i = 0; i < cloth->mConstraintRestLengths.size(); ++i)
			g_buffers->springLengths.push_back(cloth->mConstraintRestLengths[i]);

		mCloths.push_back(cloth);

		// add inflatable params
		g_buffers->inflatableVolumes.push_back(cloth->mRestVolume);
		g_buffers->inflatableCoefficients.push_back(cloth->mConstraintScale);
	}

	void Initialize()
	{
		mCloths.resize(0);

		float minSize = 0.75f;
		float maxSize = 1.0f;

		float radius = 0.12f;
		int group = 0;

		g_sceneLower = Vec3(0.0f, 0.0f, -2.5f);
		g_sceneUpper = Vec3(10.0f, 5.0f, -2.5f);

		const char* meshes[2] =
		{
			"../../data/box_high_weld.ply",
			"../../data/sphere.ply"
		};

		mPressure = 1.35f; // Overpressure

		for (int y = 1; y < 2; ++y)
		{
			for (int i = 0; i < 1; ++i)
			{
				Mesh* mesh = ImportMesh(GetFilePathByPlatform(meshes[(i + y) & 1]).c_str());
				mesh->Normalize();
				//mesh->Transform(TranslationMatrix(Point3(i*2.0f, 1.0f + y * 2.0f, 1.5f)));
				mesh->Transform(TranslationMatrix(Point3(5.0f, 1.0f + y * 2.0f, 1.5f)));

				AddCell(mesh, mPressure, NvFlexMakePhase(group++, 0));

				delete mesh;
			}
		}

		g_params.radius = radius;
		g_params.dynamicFriction = 0.4f;
		g_params.dissipation = 0.0f;
		g_params.numIterations = 10;
		g_params.particleCollisionMargin = g_params.radius*0.05f;
		g_params.drag = 0.0f;
		g_params.collisionDistance = 0.01f;
		g_params.gravity[1] = 0.0f;

		// better convergence with global relaxation factor
		g_params.relaxationMode = eNvFlexRelaxationGlobal;
		g_params.relaxationFactor = 0.25f;
		g_params.numPlanes = 6;

		g_windStrength = 0.0f;

		g_numSubsteps = 2;

		// draw options		
		g_drawPoints = false;
		g_drawSprings = 0;
		g_drawCloth = false;
	}

	virtual void DoGui() // Include a parameter varying GUI
	{
		if (imguiSlider("Over Pressure", &mPressure, 0.25f, 3.0f, 0.001f))
		{
			for (int i = 0; i < int(g_buffers->inflatablePressures.size()); ++i)
				g_buffers->inflatablePressures[i] = mPressure;
		}
	}

	virtual void Update() // Update any buffers (all guaranteed to be mapped here)
	{

		// Apply instantaneous force to particles in random section of the cell surface
		// Parameters
		float x_centre = 0, y_centre = 0, z_centre = 0;
		float azimuth, zenith, dX, dY, dZ, r;
		float angle_range = 25*(3.14159/180);
		// Random direction
		float rand_azim = Randf(-3.1415f, 3.1415f);
		float rand_zeni = Randf(0.0f, 3.1415f);
		// Random velocity
		float ins_vel = Randf(0.0f, 5.0f); 
		// Find centre of mass of the cell (average all cell positions in cell)
		for (int i = 0; i < int(NvFlexGetActiveCount(g_flex)); ++i) {
			x_centre = x_centre + g_buffers->positions[i].x;
			y_centre = y_centre + g_buffers->positions[i].y;
			z_centre = z_centre + g_buffers->positions[i].z;
		}
		x_centre = x_centre / NvFlexGetActiveCount(g_flex);
		y_centre = y_centre / NvFlexGetActiveCount(g_flex);
		z_centre = z_centre / NvFlexGetActiveCount(g_flex);

		// Apply instantaneous surface velocity to random direction in cell
		for (int i = 0; i < int(NvFlexGetActiveCount(g_flex)); ++i) {

			dX = (g_buffers->positions[i].x - x_centre);
			dY = (g_buffers->positions[i].y - y_centre);
			dZ = (g_buffers->positions[i].z - z_centre);

			r = sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));

			azimuth = atan2(dY, dX);
			zenith = acos(dZ/r);

			if (azimuth > (rand_azim - angle_range) && azimuth < (rand_azim + angle_range) && zenith > (rand_zeni - angle_range) && zenith < (rand_zeni + angle_range)) {

				g_buffers->velocities[i].x = ins_vel * cos(azimuth)*sin(zenith);
				g_buffers->velocities[i].y = ins_vel * sin(azimuth)*sin(zenith);
				g_buffers->velocities[i].z = ins_vel * cos(zenith);

			}

		}

	}

	virtual void Sync() // Send any changes to flex (all buffers guaranteed to be unmapped here)
	{
		NvFlexSetInflatables(g_flex, g_buffers->inflatableTriOffsets.buffer, g_buffers->inflatableTriCounts.buffer, g_buffers->inflatableVolumes.buffer, g_buffers->inflatablePressures.buffer, g_buffers->inflatableCoefficients.buffer, mCloths.size());
		NvFlexSetVelocities(g_flex, g_buffers->velocities.buffer, g_buffers->velocities.size());
	}

	virtual void Draw(int pass)
	{
		if (!g_drawMesh)
			return;

		int indexStart = 0;

		for (size_t i = 0; i < mCloths.size(); ++i)
		{
			DrawCloth(&g_buffers->positions[0], &g_buffers->normals[0], NULL, &g_buffers->triangles[indexStart], mCloths[i]->mTris.size(), g_buffers->positions.size(), 2, g_params.radius*0.35f);

			indexStart += mCloths[i]->mTris.size() * 3;
		}
	}

	float mPressure;

	std::vector<ClothMesh*> mCloths;
};
