#ifndef WaterMACGrid_H_
#define WaterMACGrid_H_

#include "mac_grid.h"

class WaterMACGrid : public MACGrid {
public:
	// functions to override

	void reset();
	WaterMACGrid();
	~WaterMACGrid();
	WaterMACGrid(const WaterMACGrid& orig);
	WaterMACGrid& operator=(const WaterMACGrid& orig);

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt);
	void addExternalForces(double dt);
	void project(double dt);
	void advectTemperature(double dt);
	//void advectDensity(double dt);
	void checkDivergence();

protected:
	// Setup:
	void initialize();

	// Simulation:
	void computeGravity(double dt);
	void computeBuoyancy(double dt);
	void computeVorticityConfinement(double dt);

	//Solver 
	bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	// Sets up the A matrix:
	void setUpAMatrix();

	//Drawing:
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);


protected:
	// level set signed distance
	GridDataLSet mLSet;

};


#endif