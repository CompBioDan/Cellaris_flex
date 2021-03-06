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

#include "stdafx.h"
#include <iostream>
#include <map>
#include <cassert>
#include <time.h>
#include <random>

#include "flex/core/types.h"
#include "flex/core/maths.h"
#include "flex/core/platform.h"
#include "flex/core/mesh.h"
#include "flex/core/voxelize.h"
#include "flex/core/sdf.h"
#include "flex/core/pfm.h"
#include "flex/core/tga.h"
#include "flex/core/perlin.h"
#include "flex/core/convex.h"
#include "flex/core/cloth.h"

//#include "flex/external/SDL2-2.0.4/include/SDL.h"

#include "flex/include/NvFlex.h"
#include "flex/include/NvFlexExt.h"
#include "flex/include/NvFlexDevice.h"

//#include "flex/flexSetup.h"

//#include "flex/flexHelper.h"

#include "cells/cell.h"
#include "scenes/scenes.h"
#include "utilities/scenetime.h"
#include "cells/cellpopulation.h"

#include "shaders.h"
#include "imgui.h"
#include <algorithm>
#include <Windows.h>

#include "utilities/global.h"

using namespace std;


//////// FLEX SETUPS
FluidRenderer* g_fluidRenderer; // shadersGL (struct containing fluid rendering data)
FluidRenderBuffers g_fluidRenderBuffers; // Buffer for fluid render data allows flex to write directly to vertex buffer objects
DiffuseRenderBuffers g_diffuseRenderBuffers; // Additional fluid diffusing data buffers

// GLOBAL
NvFlexSolver * g_flex; // Instance of the NvFlexSolver
NvFlexLibrary* g_flexLib; // Instance of flex library
NvFlexParams g_params;
NvFlexTimers g_timers;
int g_numDetailTimers;
NvFlexDetailTimer * g_detailTimers;
bool g_Error = false;
char g_deviceName[256];

// parameters for particles in simulations
int g_maxDiffuseParticles;
unsigned char g_maxNeighborsPerParticle;
int g_numExtraParticles;
int g_numExtraMultiplier = 1;

// mesh used for deformable object rendering
Mesh* g_mesh;
vector<int> g_meshSkinIndices;
vector<float> g_meshSkinWeights;
vector<Point3> g_meshRestPositions;
const int g_numSkinWeights = 4;

// mapping of collision mesh to render mesh
std::map<NvFlexConvexMeshId, GpuMesh*> g_convexes;
std::map<NvFlexTriangleMeshId, GpuMesh*> g_meshes;
std::map<NvFlexDistanceFieldId, GpuMesh*> g_fields;

class CreateScene;
vector<CreateScene*> g_scenes;


struct SimBuffers {

	NvFlexVector<Vec4> positions;
	NvFlexVector<Vec4> restPositions;
	NvFlexVector<Vec3> velocities;
	NvFlexVector<int> phases;
	NvFlexVector<float> densities;
	NvFlexVector<Vec4> anisotropy1;
	NvFlexVector<Vec4> anisotropy2;
	NvFlexVector<Vec4> anisotropy3;
	NvFlexVector<Vec4> normals;
	NvFlexVector<Vec4> smoothPositions;
	NvFlexVector<Vec4> diffusePositions;
	NvFlexVector<Vec4> diffuseVelocities;
	NvFlexVector<int> diffuseIndices;
	NvFlexVector<int> activeIndices;

	// convexes
	NvFlexVector<NvFlexCollisionGeometry> shapeGeometry;
	NvFlexVector<Vec4> shapePositions;
	NvFlexVector<Quat> shapeRotations;
	NvFlexVector<Vec4> shapePrevPositions;
	NvFlexVector<Quat> shapePrevRotations;
	NvFlexVector<int> shapeFlags;

	// rigids
	NvFlexVector<int> rigidOffsets;
	NvFlexVector<int> rigidIndices;
	NvFlexVector<int> rigidMeshSize;
	NvFlexVector<float> rigidCoefficients;
	NvFlexVector<Quat> rigidRotations;
	NvFlexVector<Vec3> rigidTranslations;
	NvFlexVector<Vec3> rigidLocalPositions;
	NvFlexVector<Vec4> rigidLocalNormals;

	// inflatables
	NvFlexVector<int> inflatableTriOffsets;
	NvFlexVector<int> inflatableTriCounts;
	NvFlexVector<float> inflatableVolumes;
	NvFlexVector<float> inflatableCoefficients;
	NvFlexVector<float> inflatablePressures;

	// springs
	NvFlexVector<int> springIndices;
	NvFlexVector<float> springLengths;
	NvFlexVector<float> springStiffness;

	NvFlexVector<int> triangles;
	NvFlexVector<Vec3> triangleNormals;
	NvFlexVector<Vec3> uvs;

	SimBuffers(NvFlexLibrary* l) :
		positions(l), restPositions(l), velocities(l), phases(l), densities(l),
		anisotropy1(l), anisotropy2(l), anisotropy3(l), normals(l), smoothPositions(l),
		diffusePositions(l), diffuseVelocities(l), diffuseIndices(l), activeIndices(l),
		shapeGeometry(l), shapePositions(l), shapeRotations(l), shapePrevPositions(l),
		shapePrevRotations(l), shapeFlags(l), rigidOffsets(l), rigidIndices(l), rigidMeshSize(l),
		rigidCoefficients(l), rigidRotations(l), rigidTranslations(l),
		rigidLocalPositions(l), rigidLocalNormals(l), inflatableTriOffsets(l),
		inflatableTriCounts(l), inflatableVolumes(l), inflatableCoefficients(l),
		inflatablePressures(l), springIndices(l), springLengths(l),
		springStiffness(l), triangles(l), triangleNormals(l), uvs(l)
	{}

};

SimBuffers* g_buffers;

void MapBuffers(SimBuffers* buffers)
{
	buffers->positions.map();
	buffers->restPositions.map();
	buffers->velocities.map();
	buffers->phases.map();
	buffers->densities.map();
	buffers->anisotropy1.map();
	buffers->anisotropy2.map();
	buffers->anisotropy3.map();
	buffers->normals.map();
	buffers->diffusePositions.map();
	buffers->diffuseVelocities.map();
	buffers->diffuseIndices.map();
	buffers->smoothPositions.map();
	buffers->activeIndices.map();

	// convexes
	buffers->shapeGeometry.map();
	buffers->shapePositions.map();
	buffers->shapeRotations.map();
	buffers->shapePrevPositions.map();
	buffers->shapePrevRotations.map();
	buffers->shapeFlags.map();

	buffers->rigidOffsets.map();
	buffers->rigidIndices.map();
	buffers->rigidMeshSize.map();
	buffers->rigidCoefficients.map();
	buffers->rigidRotations.map();
	buffers->rigidTranslations.map();
	buffers->rigidLocalPositions.map();
	buffers->rigidLocalNormals.map();

	buffers->springIndices.map();
	buffers->springLengths.map();
	buffers->springStiffness.map();

	// inflatables
	buffers->inflatableTriOffsets.map();
	buffers->inflatableTriCounts.map();
	buffers->inflatableVolumes.map();
	buffers->inflatableCoefficients.map();
	buffers->inflatablePressures.map();

	buffers->triangles.map();
	buffers->triangleNormals.map();
	buffers->uvs.map();
}
//
//void UnmapBuffers(SimBuffers* buffers)
//{
//	// particles
//	buffers->positions.unmap();
//	buffers->restPositions.unmap();
//	buffers->velocities.unmap();
//	buffers->phases.unmap();
//	buffers->densities.unmap();
//	buffers->anisotropy1.unmap();
//	buffers->anisotropy2.unmap();
//	buffers->anisotropy3.unmap();
//	buffers->normals.unmap();
//	buffers->diffusePositions.unmap();
//	buffers->diffuseVelocities.unmap();
//	buffers->diffuseIndices.unmap();
//	buffers->smoothPositions.unmap();
//	buffers->activeIndices.unmap();
//
//	// convexes
//	buffers->shapeGeometry.unmap();
//	buffers->shapePositions.unmap();
//	buffers->shapeRotations.unmap();
//	buffers->shapePrevPositions.unmap();
//	buffers->shapePrevRotations.unmap();
//	buffers->shapeFlags.unmap();
//
//	// rigids
//	buffers->rigidOffsets.unmap();
//	buffers->rigidIndices.unmap();
//	buffers->rigidMeshSize.unmap();
//	buffers->rigidCoefficients.unmap();
//	buffers->rigidRotations.unmap();
//	buffers->rigidTranslations.unmap();
//	buffers->rigidLocalPositions.unmap();
//	buffers->rigidLocalNormals.unmap();
//
//	// springs
//	buffers->springIndices.unmap();
//	buffers->springLengths.unmap();
//	buffers->springStiffness.unmap();
//
//	// inflatables
//	buffers->inflatableTriOffsets.unmap();
//	buffers->inflatableTriCounts.unmap();
//	buffers->inflatableVolumes.unmap();
//	buffers->inflatableCoefficients.unmap();
//	buffers->inflatablePressures.unmap();
//
//	// triangles
//	buffers->triangles.unmap();
//	buffers->triangleNormals.unmap();
//	buffers->uvs.unmap();
//}
//
//SimBuffers* AllocBuffers(NvFlexLibrary* lib)
//{
//	return new SimBuffers(lib);
//}
//
//void DestroyBuffers(SimBuffers* buffers)
//{
//	// particles
//	buffers->positions.destroy();
//	buffers->restPositions.destroy();
//	buffers->velocities.destroy();
//	buffers->phases.destroy();
//	buffers->densities.destroy();
//	buffers->anisotropy1.destroy();
//	buffers->anisotropy2.destroy();
//	buffers->anisotropy3.destroy();
//	buffers->normals.destroy();
//	buffers->diffusePositions.destroy();
//	buffers->diffuseVelocities.destroy();
//	buffers->diffuseIndices.destroy();
//	buffers->smoothPositions.destroy();
//	buffers->activeIndices.destroy();
//
//	// convexes
//	buffers->shapeGeometry.destroy();
//	buffers->shapePositions.destroy();
//	buffers->shapeRotations.destroy();
//	buffers->shapePrevPositions.destroy();
//	buffers->shapePrevRotations.destroy();
//	buffers->shapeFlags.destroy();
//
//	// rigids
//	buffers->rigidOffsets.destroy();
//	buffers->rigidIndices.destroy();
//	buffers->rigidMeshSize.destroy();
//	buffers->rigidCoefficients.destroy();
//	buffers->rigidRotations.destroy();
//	buffers->rigidTranslations.destroy();
//	buffers->rigidLocalPositions.destroy();
//	buffers->rigidLocalNormals.destroy();
//
//	// springs
//	buffers->springIndices.destroy();
//	buffers->springLengths.destroy();
//	buffers->springStiffness.destroy();
//
//	// inflatables
//	buffers->inflatableTriOffsets.destroy();
//	buffers->inflatableTriCounts.destroy();
//	buffers->inflatableVolumes.destroy();
//	buffers->inflatableCoefficients.destroy();
//	buffers->inflatablePressures.destroy();
//
//	// triangles
//	buffers->triangles.destroy();
//	buffers->triangleNormals.destroy();
//	buffers->uvs.destroy();
//
//	delete buffers;
//}
//
//struct Emitter
//{
//	Emitter() : mSpeed(0.0f), mEnabled(false), mLeftOver(0.0f), mWidth(8) {}
//
//	Vec3 mPos;
//	Vec3 mDir;
//	Vec3 mRight;
//	float mSpeed;
//	bool mEnabled;
//	float mLeftOver;
//	int mWidth;
//};
//
//vector<Emitter> g_emitters(1);	// first emitter is the camera 'gun'
//
//struct Rope
//{
//	std::vector<int> mIndices;
//};
//
//vector<Rope> g_ropes;
//
//// Viewing/camera parameters
//Vec3 g_camPos(6.0f, 8.0f, 18.0f);
//Vec3 g_camAngle(0.0f, -DegToRad(20.0f), 0.0f);
//Vec3 g_camVel(0.0f);
//Vec3 g_camSmoothVel(0.0f);
//
//float g_camSpeed;
//float g_camNear;
//float g_camFar;
//
//Vec3 g_lightPos;
//Vec3 g_lightDir;
//Vec3 g_lightTarget;


//#include "flex/include/helpers.h"
//#include "flex/include/createscenes.h"
//#include "flex/include/benchmark.h"

//void ErrorCallback(NvFlexErrorSeverity, const char* msg, const char* file, int line)
//{
//	printf("Flex: %s - %s:%d\n", msg, file, line);
//	g_Error = true;
//	//assert(0); asserts are bad for TeamCity
//}

//void CalculateRigidLocalPositions(const Vec4* restPositions, int numRestPositions, const int* offsets, const int* indices, int numRigids, Vec3* localPositions)
//{
//
//	// To improve the accuracy of the result, first transform the restPositions to relative coordinates (by finding the mean and subtracting that from all points)
//	// Note: If this is not done, one might see ghost forces if the mean of the restPositions is far from the origin.
//
//	// Calculate mean
//	Vec3 shapeOffset(0.0f);
//
//	for (int i = 0; i < numRestPositions; i++)
//	{
//		shapeOffset += Vec3(restPositions[i]);
//	}
//
//	shapeOffset /= float(numRestPositions);
//
//	int count = 0;
//
//	for (int r = 0; r < numRigids; ++r)
//	{
//		const int startIndex = offsets[r];
//		const int endIndex = offsets[r + 1];
//
//		const int n = endIndex - startIndex;
//
//		assert(n);
//
//		Vec3 com;
//
//		for (int i = startIndex; i < endIndex; ++i)
//		{
//			const int r = indices[i];
//
//			// By substracting meshOffset the calculation is done in relative coordinates
//			com += Vec3(restPositions[r]) - shapeOffset;
//		}
//
//		com /= float(n);
//
//		for (int i = startIndex; i < endIndex; ++i)
//		{
//			const int r = indices[i];
//
//			// By substracting meshOffset the calculation is done in relative coordinates
//			localPositions[count++] = (Vec3(restPositions[r]) - shapeOffset) - com;
//		}
//	}
//}

//void GetParticleBounds(Vec3& lower, Vec3& upper)
//{
//	lower = Vec3(FLT_MAX);
//	upper = Vec3(-FLT_MAX);
//
//	for (int i = 0; i < g_buffers->positions.size(); ++i)
//	{
//		lower = Min(Vec3(g_buffers->positions[i]), lower);
//		upper = Max(Vec3(g_buffers->positions[i]), upper);
//	}
//}

//void GetShapeBounds(Vec3& totalLower, Vec3& totalUpper)
//{
//	Bounds totalBounds;
//
//	for (int i = 0; i < g_buffers->shapeFlags.size(); ++i)
//	{
//		NvFlexCollisionGeometry geo = g_buffers->shapeGeometry[i];
//
//		int type = g_buffers->shapeFlags[i] & eNvFlexShapeFlagTypeMask;
//
//		Vec3 localLower;
//		Vec3 localUpper;
//
//		switch (type)
//		{
//		case eNvFlexShapeBox:
//		{
//			localLower = -Vec3(geo.box.halfExtents);
//			localUpper = Vec3(geo.box.halfExtents);
//			break;
//		}
//		case eNvFlexShapeSphere:
//		{
//			localLower = -geo.sphere.radius;
//			localUpper = geo.sphere.radius;
//			break;
//		}
//		case eNvFlexShapeCapsule:
//		{
//			localLower = -Vec3(geo.capsule.halfHeight, 0.0f, 0.0f) - Vec3(geo.capsule.radius);
//			localUpper = Vec3(geo.capsule.halfHeight, 0.0f, 0.0f) + Vec3(geo.capsule.radius);
//			break;
//		}
//		case eNvFlexShapeConvexMesh:
//		{
//			NvFlexGetConvexMeshBounds(g_flexLib, geo.convexMesh.mesh, localLower, localUpper);
//
//			// apply instance scaling
//			localLower *= geo.convexMesh.scale;
//			localUpper *= geo.convexMesh.scale;
//			break;
//		}
//		case eNvFlexShapeTriangleMesh:
//		{
//			NvFlexGetTriangleMeshBounds(g_flexLib, geo.triMesh.mesh, localLower, localUpper);
//
//			// apply instance scaling
//			localLower *= Vec3(geo.triMesh.scale);
//			localUpper *= Vec3(geo.triMesh.scale);
//			break;
//		}
//		case eNvFlexShapeSDF:
//		{
//			localLower = 0.0f;
//			localUpper = geo.sdf.scale;
//			break;
//		}
//		};
//
//		// transform local bounds to world space
//		Vec3 worldLower, worldUpper;
//		TransformBounds(localLower, localUpper, Vec3(g_buffers->shapePositions[i]), g_buffers->shapeRotations[i], 1.0f, worldLower, worldUpper);
//
//		totalBounds = Union(totalBounds, Bounds(worldLower, worldUpper));
//	}
//
//	totalLower = totalBounds.lower;
//	totalUpper = totalBounds.upper;
//}


//void Init(int scene, bool centerCamera = true)
//{
//	RandInit(); // Seed random number?
//
//	if (g_flex) // If instance of solver already exists remove it
//	{
//		if (g_buffers)
//			DestroyBuffers(g_buffers);
//
//		for (auto& iter : g_meshes)
//		{
//			NvFlexDestroyTriangleMesh(g_flexLib, iter.first);
//			DestroyGpuMesh(iter.second);
//		}
//
//		for (auto& iter : g_fields)
//		{
//			NvFlexDestroyDistanceField(g_flexLib, iter.first);
//			DestroyGpuMesh(iter.second);
//		}
//
//		for (auto& iter : g_convexes)
//		{
//			NvFlexDestroyConvexMesh(g_flexLib, iter.first);
//			DestroyGpuMesh(iter.second);
//		}
//
//
//		g_fields.clear();
//		g_meshes.clear();
//		g_convexes.clear();
//
//		NvFlexDestroySolver(g_flex);
//		g_flex = NULL;
//	}
//
//	// alloc buffers
//	g_buffers = AllocBuffers(g_flexLib);
//
//	// map during initialization
//	MapBuffers(g_buffers);
//
//	g_buffers->positions.resize(0);
//	g_buffers->velocities.resize(0);
//	g_buffers->phases.resize(0);
//
//	g_buffers->rigidOffsets.resize(0);
//	g_buffers->rigidIndices.resize(0);
//	g_buffers->rigidMeshSize.resize(0);
//	g_buffers->rigidRotations.resize(0);
//	g_buffers->rigidTranslations.resize(0);
//	g_buffers->rigidCoefficients.resize(0);
//	g_buffers->rigidLocalPositions.resize(0);
//	g_buffers->rigidLocalNormals.resize(0);
//
//	g_buffers->springIndices.resize(0);
//	g_buffers->springLengths.resize(0);
//	g_buffers->springStiffness.resize(0);
//	g_buffers->triangles.resize(0);
//	g_buffers->triangleNormals.resize(0);
//	g_buffers->uvs.resize(0);
//
//	g_meshSkinIndices.resize(0);
//	g_meshSkinWeights.resize(0);
//
//	g_emitters.resize(1);
//	g_emitters[0].mEnabled = false;
//	g_emitters[0].mSpeed = 1.0f;
//
//	g_buffers->shapeGeometry.resize(0);
//	g_buffers->shapePositions.resize(0);
//	g_buffers->shapeRotations.resize(0);
//	g_buffers->shapePrevPositions.resize(0);
//	g_buffers->shapePrevRotations.resize(0);
//	g_buffers->shapeFlags.resize(0);
//
//	g_ropes.resize(0);
//
//	// remove collision shapes
//	delete g_mesh; g_mesh = NULL;
//
//	g_frame = 0;
//	g_pause = false;
//
//	g_dt = 1.0f / 60.0f;
//	g_waveTime = 0.0f;
//	g_windTime = 0.0f;
//	g_windStrength = 1.0f;
//
//	g_blur = 1.0f;
//	XVector4<float> g_fluidColor = Vec4(0.1f, 0.4f, 0.8f, 1.0f);
//	XVector3<float> g_meshColor = Vec3(0.9f, 0.9f, 0.9f);
//	g_drawEllipsoids = false;
//	g_drawPoints = true;
//	g_drawCloth = true;
//	g_expandCloth = 0.0f;
//
//	g_drawOpaque = false;
//	g_drawSprings = false;
//	g_drawDiffuse = false;
//	g_drawMesh = true;
//	g_drawRopes = true;
//	g_drawDensity = false;
//	g_ior = 1.0f;
//	g_lightDistance = 2.0f;
//	g_fogDistance = 0.005f;
//
//	g_camSpeed = 0.075f;
//	g_camNear = 0.01f;
//	g_camFar = 1000.0f;
//
//	g_pointScale = 1.0f;
//	g_ropeScale = 1.0f;
//	g_drawPlaneBias = 0.0f;
//
//	// sim params
//	g_params.gravity[0] = 0.0f;
//	g_params.gravity[1] = -9.8f;
//	g_params.gravity[2] = 0.0f;
//
//	g_params.wind[0] = 0.0f;
//	g_params.wind[1] = 0.0f;
//	g_params.wind[2] = 0.0f;
//
//	g_params.radius = 0.15f;
//	g_params.viscosity = 0.0f;
//	g_params.dynamicFriction = 0.0f;
//	g_params.staticFriction = 0.0f;
//	g_params.particleFriction = 0.0f; // scale friction between particles by default
//	g_params.freeSurfaceDrag = 0.0f;
//	g_params.drag = 0.0f;
//	g_params.lift = 0.0f;
//	g_params.numIterations = 3;
//	g_params.fluidRestDistance = 0.0f;
//	g_params.solidRestDistance = 0.0f;
//
//	g_params.anisotropyScale = 1.0f;
//	g_params.anisotropyMin = 0.1f;
//	g_params.anisotropyMax = 2.0f;
//	g_params.smoothing = 1.0f;
//
//	g_params.dissipation = 0.0f;
//	g_params.damping = 0.0f;
//	g_params.particleCollisionMargin = 0.0f;
//	g_params.shapeCollisionMargin = 0.0f;
//	g_params.collisionDistance = 0.0f;
//	g_params.plasticThreshold = 0.0f;
//	g_params.plasticCreep = 0.0f;
//	g_params.fluid = false;
//	g_params.sleepThreshold = 0.0f;
//	g_params.shockPropagation = 0.0f;
//	g_params.restitution = 0.0f;
//
//	g_params.maxSpeed = FLT_MAX;
//	g_params.maxAcceleration = 100.0f;	// approximately 10x gravity
//
//	g_params.relaxationMode = eNvFlexRelaxationLocal;
//	g_params.relaxationFactor = 1.0f;
//	g_params.solidPressure = 1.0f;
//	g_params.adhesion = 0.0f;
//	g_params.cohesion = 0.025f;
//	g_params.surfaceTension = 0.0f;
//	g_params.vorticityConfinement = 0.0f;
//	g_params.buoyancy = 1.0f;
//	g_params.diffuseThreshold = 100.0f;
//	g_params.diffuseBuoyancy = 1.0f;
//	g_params.diffuseDrag = 0.8f;
//	g_params.diffuseBallistic = 16;
//	g_params.diffuseSortAxis[0] = 0.0f;
//	g_params.diffuseSortAxis[1] = 0.0f;
//	g_params.diffuseSortAxis[2] = 0.0f;
//	g_params.diffuseLifetime = 2.0f;
//
//	g_numSubsteps = 2;
//
//	// planes created after particles
//	g_params.numPlanes = 1;
//
//	g_diffuseScale = 0.5f;
//	XVector4<float> g_diffuseColor = 1.0f;
//	g_diffuseMotionScale = 1.0f;
//	g_diffuseShadow = false;
//	g_diffuseInscatter = 0.8f;
//	g_diffuseOutscatter = 0.53f;
//
//	// reset phase 0 particle color to blue
//	extern Colour gColors[];
//	gColors[0] = Colour(0.0f, 0.5f, 1.0f);
//
//	g_numSolidParticles = 0;
//
//	g_waveFrequency = 1.5f;
//	g_waveAmplitude = 1.5f;
//	g_waveFloorTilt = 0.0f;
//	g_emit = false;
//	g_warmup = false;
//
//	g_mouseParticle = -1;
//
//	g_maxDiffuseParticles = 0;	// number of diffuse particles
//	g_maxNeighborsPerParticle = 96;
//	g_numExtraParticles = 0;	// number of particles allocated but not made active	
//
//	XVector3<float> g_sceneLower = FLT_MAX;
//	XVector3<float> g_sceneUpper = -FLT_MAX;
//
//	// create scene
//	g_scenes[g_scene]->Initialize();
//
//	uint32_t numParticles = g_buffers->positions.size();
//	uint32_t maxParticles = numParticles + g_numExtraParticles * g_numExtraMultiplier;
//
//	// by default solid particles use the maximum radius
//	if (g_params.fluid && g_params.solidRestDistance == 0.0f)
//		g_params.solidRestDistance = g_params.fluidRestDistance;
//	else
//		g_params.solidRestDistance = g_params.radius;
//
//	// collision distance with shapes half the radius
//	if (g_params.collisionDistance == 0.0f)
//	{
//		g_params.collisionDistance = g_params.radius*0.5f;
//
//		if (g_params.fluid)
//			g_params.collisionDistance = g_params.fluidRestDistance*0.5f;
//	}
//
//	// default particle friction to 10% of shape friction
//	if (g_params.particleFriction == 0.0f)
//		g_params.particleFriction = g_params.dynamicFriction*0.1f;
//
//	// add a margin for detecting contacts between particles and shapes
//	if (g_params.shapeCollisionMargin == 0.0f)
//		g_params.shapeCollisionMargin = g_params.collisionDistance*0.5f;
//
//	// calculate particle bounds
//	Vec3 particleLower, particleUpper;
//	GetParticleBounds(particleLower, particleUpper);
//
//	// accommodate shapes
//	Vec3 shapeLower, shapeUpper;
//	GetShapeBounds(shapeLower, shapeUpper);
//
//	// update bounds
//	g_sceneLower = Min(Min(g_sceneLower, particleLower), shapeLower);
//	g_sceneUpper = Max(Max(g_sceneUpper, particleUpper), shapeUpper);
//
//	g_sceneLower -= g_params.collisionDistance;
//	g_sceneUpper += g_params.collisionDistance;
//
//	// update collision planes to match flexs
//	Vec3 up = Normalize(Vec3(-g_waveFloorTilt, 1.0f, 0.0f));
//
//	(Vec4&)g_params.planes[0] = Vec4(up.x, up.y, up.z, 0.0f);
//	(Vec4&)g_params.planes[1] = Vec4(0.0f, 0.0f, 1.0f, -g_sceneLower.z);
//	(Vec4&)g_params.planes[2] = Vec4(1.0f, 0.0f, 0.0f, -g_sceneLower.x);
//	(Vec4&)g_params.planes[3] = Vec4(-1.0f, 0.0f, 0.0f, g_sceneUpper.x);
//	(Vec4&)g_params.planes[4] = Vec4(0.0f, 0.0f, -1.0f, g_sceneUpper.z);
//	(Vec4&)g_params.planes[5] = Vec4(0.0f, -1.0f, 0.0f, g_sceneUpper.y);
//
//	float g_wavePlane = g_params.planes[2][3];
//
//	g_buffers->diffusePositions.resize(g_maxDiffuseParticles);
//	g_buffers->diffuseVelocities.resize(g_maxDiffuseParticles);
//	g_buffers->diffuseIndices.resize(g_maxDiffuseParticles);
//
//	// for fluid rendering these are the Laplacian smoothed positions
//	g_buffers->smoothPositions.resize(maxParticles);
//
//	g_buffers->normals.resize(0);
//	g_buffers->normals.resize(maxParticles);
//
//	// initialize normals (just for rendering before simulation starts)
//	int numTris = g_buffers->triangles.size() / 3;
//	for (int i = 0; i < numTris; ++i)
//	{
//		Vec3 v0 = Vec3(g_buffers->positions[g_buffers->triangles[i * 3 + 0]]);
//		Vec3 v1 = Vec3(g_buffers->positions[g_buffers->triangles[i * 3 + 1]]);
//		Vec3 v2 = Vec3(g_buffers->positions[g_buffers->triangles[i * 3 + 2]]);
//
//		Vec3 n = Cross(v1 - v0, v2 - v0);
//
//		g_buffers->normals[g_buffers->triangles[i * 3 + 0]] += Vec4(n, 0.0f);
//		g_buffers->normals[g_buffers->triangles[i * 3 + 1]] += Vec4(n, 0.0f);
//		g_buffers->normals[g_buffers->triangles[i * 3 + 2]] += Vec4(n, 0.0f);
//	}
//
//	for (int i = 0; i < int(maxParticles); ++i)
//		g_buffers->normals[i] = Vec4(SafeNormalize(Vec3(g_buffers->normals[i]), Vec3(0.0f, 1.0f, 0.0f)), 0.0f);
//
//
//	// save mesh positions for skinning
//	if (g_mesh)
//	{
//		g_meshRestPositions = g_mesh->m_positions;
//	}
//	else
//	{
//		g_meshRestPositions.resize(0);
//	}
//
//	// main create method for the Flex solver
//	g_flex = NvFlexCreateSolver(g_flexLib, maxParticles, g_maxDiffuseParticles, g_maxNeighborsPerParticle);
//
//	// give scene a chance to do some post solver initialization
//	g_scenes[g_scene]->PostInitialize();
//
//	// center camera on particles
//	if (centerCamera)
//	{
//		g_camPos = Vec3((g_sceneLower.x + g_sceneUpper.x)*0.5f, min(g_sceneUpper.y*1.25f, 6.0f), g_sceneUpper.z + min(g_sceneUpper.y, 6.0f)*2.0f);
//		g_camAngle = Vec3(0.0f, -DegToRad(15.0f), 0.0f);
//
//		// give scene a chance to modify camera position
//		g_scenes[g_scene]->CenterCamera();
//	}
//
//	// create active indices (just a contiguous block for the demo)
//	g_buffers->activeIndices.resize(g_buffers->positions.size());
//	for (int i = 0; i < g_buffers->activeIndices.size(); ++i)
//		g_buffers->activeIndices[i] = i;
//
//	// resize particle buffers to fit
//	g_buffers->positions.resize(maxParticles);
//	g_buffers->velocities.resize(maxParticles);
//	g_buffers->phases.resize(maxParticles);
//
//	g_buffers->densities.resize(maxParticles);
//	g_buffers->anisotropy1.resize(maxParticles);
//	g_buffers->anisotropy2.resize(maxParticles);
//	g_buffers->anisotropy3.resize(maxParticles);
//
//	// save rest positions
//	g_buffers->restPositions.resize(g_buffers->positions.size());
//	for (int i = 0; i < g_buffers->positions.size(); ++i)
//		g_buffers->restPositions[i] = g_buffers->positions[i];
//
//	// builds rigids constraints
//	if (g_buffers->rigidOffsets.size())
//	{
//		assert(g_buffers->rigidOffsets.size() > 1);
//
//		const int numRigids = g_buffers->rigidOffsets.size() - 1;
//
//		// calculate local rest space positions
//		g_buffers->rigidLocalPositions.resize(g_buffers->rigidOffsets.back());
//		CalculateRigidLocalPositions(&g_buffers->positions[0], g_buffers->positions.size(), &g_buffers->rigidOffsets[0], &g_buffers->rigidIndices[0], numRigids, &g_buffers->rigidLocalPositions[0]);
//
//		g_buffers->rigidRotations.resize(g_buffers->rigidOffsets.size() - 1, Quat());
//		g_buffers->rigidTranslations.resize(g_buffers->rigidOffsets.size() - 1, Vec3());
//
//	}
//
//	// unmap so we can start transferring data to GPU
//	UnmapBuffers(g_buffers);
//
//	//-----------------------------
//	// Send data to Flex
//
//	NvFlexSetParams(g_flex, &g_params);
//	NvFlexSetParticles(g_flex, g_buffers->positions.buffer, numParticles);
//	NvFlexSetVelocities(g_flex, g_buffers->velocities.buffer, numParticles);
//	NvFlexSetNormals(g_flex, g_buffers->normals.buffer, numParticles);
//	NvFlexSetPhases(g_flex, g_buffers->phases.buffer, g_buffers->phases.size());
//	NvFlexSetRestParticles(g_flex, g_buffers->restPositions.buffer, g_buffers->restPositions.size());
//
//	NvFlexSetActive(g_flex, g_buffers->activeIndices.buffer, numParticles);
//
//	// springs
//	if (g_buffers->springIndices.size())
//	{
//		assert((g_buffers->springIndices.size() & 1) == 0);
//		assert((g_buffers->springIndices.size() / 2) == g_buffers->springLengths.size());
//
//		NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
//	}
//
//	// rigids
//	if (g_buffers->rigidOffsets.size())
//	{
//		NvFlexSetRigids(g_flex, g_buffers->rigidOffsets.buffer, g_buffers->rigidIndices.buffer, g_buffers->rigidLocalPositions.buffer, g_buffers->rigidLocalNormals.buffer, g_buffers->rigidCoefficients.buffer, g_buffers->rigidRotations.buffer, g_buffers->rigidTranslations.buffer, g_buffers->rigidOffsets.size() - 1, g_buffers->rigidIndices.size());
//	}
//
//	// inflatables
//	if (g_buffers->inflatableTriOffsets.size())
//	{
//		NvFlexSetInflatables(g_flex, g_buffers->inflatableTriOffsets.buffer, g_buffers->inflatableTriCounts.buffer, g_buffers->inflatableVolumes.buffer, g_buffers->inflatablePressures.buffer, g_buffers->inflatableCoefficients.buffer, g_buffers->inflatableTriOffsets.size());
//	}
//
//	// dynamic triangles
//	if (g_buffers->triangles.size())
//	{
//		NvFlexSetDynamicTriangles(g_flex, g_buffers->triangles.buffer, g_buffers->triangleNormals.buffer, g_buffers->triangles.size() / 3);
//	}
//
//	// collision shapes
//	if (g_buffers->shapeFlags.size())
//	{
//		NvFlexSetShapes(
//			g_flex,
//			g_buffers->shapeGeometry.buffer,
//			g_buffers->shapePositions.buffer,
//			g_buffers->shapeRotations.buffer,
//			g_buffers->shapePrevPositions.buffer,
//			g_buffers->shapePrevRotations.buffer,
//			g_buffers->shapeFlags.buffer,
//			int(g_buffers->shapeFlags.size()));
//	}
//
//	// create render buffers
//	g_fluidRenderBuffers = CreateFluidRenderBuffers(maxParticles, g_interop);
//	g_diffuseRenderBuffers = CreateDiffuseRenderBuffers(g_maxDiffuseParticles, g_interop);
//
//	// perform initial sim warm up
//	if (g_warmup)
//	{
//		printf("Warming up sim..\n");
//
//		// warm it up (relax positions to reach rest density without affecting velocity)
//		NvFlexParams copy = g_params;
//		copy.numIterations = 4;
//
//		NvFlexSetParams(g_flex, &copy);
//
//		const int kWarmupIterations = 100;
//
//		for (int i = 0; i < kWarmupIterations; ++i)
//		{
//			NvFlexUpdateSolver(g_flex, 0.0001f, 1, false);
//			NvFlexSetVelocities(g_flex, g_buffers->velocities.buffer, maxParticles);
//		}
//
//		// udpate host copy
//		NvFlexGetParticles(g_flex, g_buffers->positions.buffer, g_buffers->positions.size());
//		NvFlexGetSmoothParticles(g_flex, g_buffers->smoothPositions.buffer, g_buffers->smoothPositions.size());
//		NvFlexGetAnisotropy(g_flex, g_buffers->anisotropy1.buffer, g_buffers->anisotropy2.buffer, g_buffers->anisotropy3.buffer);
//
//		printf("Finished warm up.\n");
//	}
//}



int main()
{

	//// Start of the testing for the Cellaris simulating framework
	std::cout << "Starting testing for the Cellaris framework..." << '\n' << '\n';

	// seed random number for random age for initial cells
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> dis(-8, -4);
	std::uniform_real_distribution<double> pos(0, 4);

	// Set up the 'scene_time' instance (containing the start time, end time and timesteps)
	double end_time = 25.0; 
	int time_steps = 250;
	int num_cells = 1;

	// Create new scene time instance and set initial time, end time and timesteps
	SceneTime* p_scene_time = SceneTime::Instance();
	p_scene_time->SetStartTime(0.0);
	p_scene_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);

	// Generate a population of cells (stored within a vector of Cell*)
	std::vector<Cell*> cell_population; // vector containing the cell information
	cell_population.reserve(num_cells); // reserve enough space for Cell*s

	// CHANGE: initialising positions for the cells, eventually change this to allow for setting cell positions or picking random within domain
	doubleVec3d nPosition; 

	int particleoffset = 0;

	for (int i = 0; i < num_cells; i++)
	{
		Cell* p_cell(new Cell());

		// Allocate a birth time for the cell
		p_cell->setBirthTime(dis(gen));

		// Allocate an initial placement of the cell
		nPosition.x = (1.0 + pos(gen)); nPosition.y = (1.0 + pos(gen)); nPosition.z = (1.0 + pos(gen));
		p_cell->setCellPos(nPosition);

		// Allocate cell id
		p_cell->setCellId(i);

		// Allocate number of particles in the cell (CHANGE: in future this will depend on specific cell-type)
		p_cell->setNumberParticles(3);

		// Needed for flex, particle offset for the buffers
		p_cell->setParticleOffset(particleoffset);

		// Add cell to population
		cell_population.push_back(p_cell);

		particleoffset += 3;
	}


	// Create our scene instance and set initial birth count to zero
	Scene* p_scene = Scene::Instance();
	p_scene->setBirths(0); // number of births at start of the scene is zero
	p_scene->setEndTime(end_time); // set end time for the scene
	p_scene->setDt(end_time / time_steps); // set the timestep
	p_scene->setNumberActiveParticles(particleoffset); // number of active particles in the scene 

	// Add the cell population to the scene
	for (unsigned i = 0; i < cell_population.size(); i++)
	{
		p_scene->addCell(cell_population[i]);
	}

	std::cout << "Number of cells in scene " << p_scene->getNumberCells() << '\n';

	p_scene->Solve();

	std::cout << "Scene solved! " << '\n';

	std::cout << "Number of cells in scene " << p_scene->getNumberCells() << '\n';

	std::cout << "number active: " << p_scene->getNumberActiveParticles() << '\n';

	std::cout << "particle offset last cell: " << p_scene->getCell(1)->getParticleOffset() << '\n';

	//// use the PhysX GPU selected from the NVIDIA control panel
	//int g_device = -1;
	//bool g_extensions = true;


	//if (g_device == -1)
	//	g_device = NvFlexDeviceGetSuggestedOrdinal();

	// Create an optimized CUDA context for Flex and set it on the 
	// calling thread. This is an optional call, it is fine to use 
	// a regular CUDA context, although creating one through this API
	// is recommended for best performance.
	//bool success = NvFlexDeviceCreateCudaContext(g_device);

	/*if (!success)
	{
		printf("Error creating CUDA context.\n");
		exit(-1);
	}*/

	//NvFlexInitDesc desc;
	//desc.deviceIndex = g_device;
	//desc.enableExtensions = g_extensions;
	//desc.renderDevice = 0;
	//desc.renderContext = 0;
	//desc.computeType = eNvFlexCUDA;

	//// Init Flex library, note that no CUDA methods should be called before this 
	//// point to ensure we get the device context we want
	//g_flexLib = NvFlexInit(NV_FLEX_VERSION, ErrorCallback, &desc);

	//if (g_Error || g_flexLib == NULL)
	//{
	//	printf("Could not initialize Flex, exiting.\n");
	//	exit(-1);
	//}

	//// store device name
	//strcpy(g_deviceName, NvFlexGetDeviceName(g_flexLib));
	//printf("Compute Device: %s\n\n", g_deviceName);

	////// init default scene
	////Init(g_scene);

	////SDLMainLoop();

	////Shutdown();

	SceneTime::Destroy();
	Scene::Destroy();
	std::cin.ignore();
	return 0;
}





