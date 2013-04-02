#ifndef WaterMACGrid_H_
#define WaterMACGrid_H_

#include "mac_grid.h"
#include <queue>

class WaterMACGrid : public MACGrid {
public:
	// display mode for signed distance
	static bool theDisplaySignedDistance;

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

	//Advect the signed distances
	void advectSignedDistances(double dt);

protected:
	// Setup:
	void initialize();

	// Simulation:
	void computeGravity(double dt);
	void computeBuoyancy(double dt);
	void computeVorticityConfinement(double dt);
	void reinitializeLevelSet();

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
	void drawSignedDistance();
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData for level set signed distance
	GridDataLSet mLSet;

};


class SignedDistCell {
public:
	SignedDistCell();
	SignedDistCell(int,int,int,double);
	~SignedDistCell();
	bool operator<(SignedDistCell &cell);

	int I;
	int J;
	int K;
	double signedDist;
};


struct GreaterDistance : public std::binary_function<SignedDistCell*, SignedDistCell*, bool>
{
    bool operator()(const SignedDistCell* first, const SignedDistCell* second) const
    {
        return first->signedDist < second->signedDist;
    }
};

#endif