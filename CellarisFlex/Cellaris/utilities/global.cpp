#include "global.h"

extern bool g_pause = false;
extern bool g_step = false;
extern bool g_capture = false;
extern bool g_showHelp = true;
extern bool g_tweakPanel = true;
extern bool g_fullscreen = false;
extern bool g_wireframe = false;
extern bool g_debug = false;

extern bool g_emit = false;
extern bool g_warmup = false;

float g_windTime = 0.0f;
float g_windFrequency = 0.1f;
float g_windStrength = 0.0f;

bool g_wavePool = false;
float g_waveTime = 0.0f;
float g_waveFrequency = 1.5f;
float g_waveAmplitude = 1.0f;
float g_waveFloorTilt = 0.0f;

bool g_drawBases = false;
bool g_drawContacts = false;
bool g_drawNormals = false;
bool g_drawShapeGrid = false;
bool g_drawDensity = false;

float g_dt = 1.0f / 60.0f;

int g_scene = 0;
int g_selectedScene = g_scene;

bool g_resetScene = false;

int g_frame = 0;
int g_numSolidParticles = 0;
int g_numSubsteps = 2;
int g_mouseParticle = -1;

bool g_interop = true;

bool g_diffuseShadow = false;
bool g_drawCloth = false;
bool g_drawDiffuse = false;
bool g_drawEllipsoids = false;
bool g_drawOpaque = false;
bool g_drawPoints = false;
bool g_drawRopes = false;
bool g_drawMesh = false;
int g_drawSprings = false;

float g_blur = 1.0f;
float g_diffuseInscatter = 0.8f;
float g_diffuseOutscatter = 0.53f;
float g_diffuseScale = 0.5f;
float g_diffuseMotionScale = 1.0f;

float g_expandCloth = 0.0f;
float g_drawPlaneBias = 0.0f;
float g_fogDistance = 0.005f;
float g_ior = 1.0f;
float g_lightDistance = 2.0f;
float g_pointScale = 1.0f;
float g_ropeScale = 1.0f;

