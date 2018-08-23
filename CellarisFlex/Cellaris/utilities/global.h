
#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "../flex/core/maths.h"

// Simulation control/parameter changing/settings
extern bool g_pause;
extern bool g_step;
extern bool g_capture;
extern bool g_showHelp;
extern bool g_tweakPanel;
extern bool g_fullscreen;
extern bool g_wireframe;
extern bool g_debug;

extern bool g_emit;
extern bool g_warmup;

extern float g_windTime;
extern float g_windFrequency;
extern float g_windStrength;

extern bool g_wavePool;
extern float g_waveTime;
extern float g_wavePlane;
extern float g_waveFrequency;
extern float g_waveAmplitude;
extern float g_waveFloorTilt;

extern Vec3 g_sceneLower;
extern Vec3 g_sceneUpper;

extern float g_blur;
extern float g_ior;
extern bool g_drawEllipsoids;
extern bool g_drawPoints;
extern bool g_drawMesh;
extern bool g_drawCloth;
extern float g_expandCloth;	// amount to expand cloth along normal (to account for particle radius)

extern bool g_drawOpaque;
extern int g_drawSprings;		// 0: no draw, 1: draw stretch 2: draw tether
extern bool g_drawBases;
extern bool g_drawContacts;
extern bool g_drawNormals;
extern bool g_drawDiffuse;
extern bool g_drawShapeGrid;
extern bool g_drawDensity;
extern bool g_drawRopes;
extern float g_pointScale;
extern float g_ropeScale;
extern float g_drawPlaneBias;	// move planes along their normal for rendering

extern float g_diffuseScale;
extern float g_diffuseMotionScale;
extern bool g_diffuseShadow;
extern float g_diffuseInscatter;
extern float g_diffuseOutscatter;

// Simulation dt
extern float g_dt;	// the time delta used for simulation
extern float g_realdt;				// the real world time delta between updates

extern float g_waitTime;		// the CPU time spent waiting for the GPU
extern float g_updateTime;     // the CPU time spent on Flex
extern float g_renderTime;		// the CPU time spent calling OpenGL to render the scene
						// the above times don't include waiting for vsync
extern float g_simLatency;     // the time the GPU spent between the first and last NvFlexUpdateSolver() operation. Because some GPUs context switch, this can include graphics time.

extern int g_scene;
extern int g_selectedScene;
extern int g_levelScroll;			// offset for level selection scroll area
extern bool g_resetScene;  //if the user clicks the reset button or presses the reset key this is set to true;

extern int g_frame;
extern int g_numSolidParticles;

extern Vec4 g_fluidColor;
extern Vec4 g_diffuseColor;
extern Vec3 g_meshColor;
extern Vec3  g_clearColor;
extern float g_lightDistance;
extern float g_fogDistance;

extern int g_numSubsteps;

extern int g_mouseParticle;

extern bool g_interop;




#endif /* GLOBAL_H_ */