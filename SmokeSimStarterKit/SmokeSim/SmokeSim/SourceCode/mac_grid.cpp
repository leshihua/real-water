// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>


// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X \
   for(int k = 0; k < theDim[MACGrid::Z]; k++) \
		for(int j = 0; j < theDim[MACGrid::Y]; j++) \
			for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_Y \
   for(int k = 0; k < theDim[MACGrid::Z]; k++) \
		for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
			for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
		for(int j = 0; j < theDim[MACGrid::Y]; j++) \
			for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define PRINT_VELOCITY \
	FOR_EACH_FACE { cout<<getVelocity(vec3(i,j,k))<<endl;}


MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT; 
   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);
   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.
	    //Top Right Corner
		mU(23,24,0) = -3.0;
		mU(24,25,0) = -3.0;
		mV(23,24,0) = -3.0;
		mV(24,25,0) = -3.0;

		//Bottom Mid Densities & Temp
		mD(17,0,0) = 1.5;
		mT(17,0,0) = 350;

		mD(18,0,0) = 1.5;
		mT(18,0,0) = 350;

		mD(19,0,0) = 1.5;
		mT(19,0,0) = 350;

		//Mid Right  
		mU(23,14,0) = -10.0;
		mU(24,15,0) = -10.0;

		mV(24,20,0) = -20;
		mV(24,19,0) = -20;
		mV(24, 18,0)=-20;

		mD(24,20,0) = 1.1;
	
		//Mid Left
		mU(20,11,0) = -3.0;
		mU(21,11,0) = -3.0;
		mV(20,11,0) = -1.0;
		mV(21,11,0) = -1.0;

		//Lower Left Side
		mD(2,10,0) = 1.5;
		mU(2,10,0) = 10;

		//Mid-Mid
		mU(10,18,0) = 10.2;
		mU(11,19,0) = 10.1;
		mU(12,20,0) = 10.2;
}

void MACGrid::advectVelocity(double dt)
{	// TODO: Calculate new velocities and store in target.
	target.mU = mU;
    target.mV = mV;
    target.mW = mW;

	//NEW WAY
	double halfCell = theCellSize/2;
	FOR_EACH_FACE {
		//OLD - WORKING WAY
		vec3 posG = getCenter(i,j,k);

		////World coordinates of X, Y, Z commponents velocities
		vec3 posGx = vec3(i*theCellSize, j*theCellSize+halfCell, k*theCellSize+halfCell);
		vec3 posGy = vec3(i*theCellSize+halfCell, j*theCellSize, k*theCellSize+halfCell);
		vec3 posGz = vec3(i*theCellSize+halfCell, j*theCellSize+halfCell, k*theCellSize);
		//vec3 xG = vec3(i,j-,k-0.5*theCellSize); 

		vec3 u = getVelocity(posGx);
		vec3 v = getVelocity(posGy);
		vec3 w = getVelocity(posGz);

		vec3 posPx = posGx-dt*u;
		vec3 posPy = posGy-dt*v;
		vec3 posPz = posGz-dt*v;

		//Get interpolated values from the grid
		vec3 tracedVel = vec3(getVelocityX(posPx), getVelocityY(posPy), getVelocityZ(posPz));	

		//Update values in the target grid
		target.mU(i,j,k) = tracedVel[0];
		target.mV(i,j,k) = tracedVel[1];
		target.mW(i,j,k) = tracedVel[2];
		
	}
    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
	//PRINT_VELOCITY
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	target.mT = mT;
    //Use RK-2 integration to figure out the advected temp at each cell pos
	FOR_EACH_CELL {
		//Current position
		vec3 xG = getCenter(i, j, k);
		//Velocity at position
		vec3 vel = getVelocity(xG);
		//Previous Position
		vec3 xP = xG-dt*vel;

		double temp = getTemperature(xP);
		target.mT(i, j, k) = temp;

	}
	
	// Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;

	FOR_EACH_CELL {
		//Current position
		vec3 xG = getCenter(i, j, k);
		//Velocity at position
		vec3 vel = getVelocity(xG);   
		//Previous Position
		vec3 xP = xG-dt*vel;
		//Update density
		double density = getDensity(xP);
		target.mD(i, j, k) = density;
	}
    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::computeBuoyancy(double dt)
{

	// TODO: Calculate bouyancy and store in target.

	double halfCell = (theCellSize/2);
	double ambientT = 300;
	double alpha = 0.11;// 0.1;
	double beta = 0.4; //0.7;

	target.mV = mV;
	//Perform bouyancy
	FOR_EACH_CELL {
		if (j > 0 && j < theDim[MACGrid::Y]) {
			vec3 currPos = getCenter(i,j,k);
			vec3 nextPos = getCenter(i,j-1,k);
			double currT, nextT, currS, nextS;	
			//Avg Temp b/w cells
			currT = mT(i,j,k);//getTemperature(currPos); //
			nextT = mT(i,j-1,k);//getTemperature(nextPos); //
			
			//Avg Concentration b/w cells
			currS = mD(i,j,k); //getDensity(currPos); //mD(i,j,k);
			nextS = mD(i,j-1,k);//getDensity(nextPos); //mD(i,j+1,k);
			
			double avgT = (currT + nextT)/2;
			double avgS = (currS + nextS)/2;

			//Apply buoyancy force to Y-component of velocity
			double FbY = (-alpha*avgS) + (beta*(avgT-ambientT));
			//vec3 FBuoy = vec3(0, FbY, 0);

			//Use Force, FbY, to get y-velocity comp. and add to this y-velocity 
			target.mV(i,j,k) = mV(i,j,k) + dt*FbY;
		}
	}
   // Then save the result to our object.
   mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{

	//Constants
	double deltaX = theCellSize;
	double dblCellSize = 2*theCellSize;		//Double cell sie
	double ep = 9.8;					//Vorticity Force Constant
	
	//Store Previous Velocities
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	//Storage for Vorticity Components & Magnitude
	GridData VortX, VortY, VortZ, VortMag;
	VortX.initialize();
	VortY.initialize(); 
	VortZ.initialize(); 
	VortMag.initialize();

	//Compute Vorticity Components & Magnitudes
	FOR_EACH_FACE {
		//Accessors for adjacent cubes
		vec3 iPrev = getCenter(i-1,j,k);
		vec3 iNext = getCenter(i+1,j,k);
				
		vec3 jPrev = getCenter(i,j-1,k);
		vec3 jNext = getCenter(i,j+1,k);
				
		vec3 kPrev = getCenter(i,j,k-1);
		vec3 kNext = getCenter(i,j,k+1);
				
		//X-gradient component
		double x1 = getVelocityZ(jNext) - getVelocityZ(jPrev);
		double x2 = getVelocityY(kNext) - getVelocityY(kPrev);
		double xVort = (x1 - x2)/dblCellSize;

		//Y-gradient component
		double y1 = getVelocityX(kNext) - getVelocityX(kPrev);
		double y2 = getVelocityZ(iNext) - getVelocityZ(iPrev);
		double yVort = (y1 - y2)/dblCellSize;

		//Z-gradient component
		double z1 = getVelocityY(iNext) - getVelocityY(iPrev);
		double z2 = getVelocityX(jNext) - getVelocityX(jPrev);
		double zVort = (z1 - z2)/dblCellSize;

		//Vorticity Vector
		vec3 vortOmega = vec3(xVort, yVort, zVort);

		//Vorticity Magnitude
		double vortMag = vortOmega.Length();

		//Store the gradient components
		VortX(i,j,k) = xVort;
		VortY(i,j,k) = yVort;
		VortZ(i,j,k) = zVort;
		//Compute VortMag
		VortMag(i,j,k) = vortMag;		
	}
	//Compute Vorticity Gradient, Normal and then Confinement Force 
	FOR_EACH_FACE {
		//Ensure this is not a boundary
		bool isNotBoundary = !(i == 0 || j == 0 || k == 0 || (i == theDim[MACGrid::X]-1) || (j == theDim[MACGrid::Y]-1) || (k == theDim[MACGrid::Z]-1));
		if (isNotBoundary) {
			double gX = (VortMag(i+1,j,k)-VortMag(i-1,j,k))/(2*theCellSize);
			double gY = (VortMag(i,j+1,k)-VortMag(i,j-1,k))/(2*theCellSize);
			double gZ = (VortMag(i,j,k+1)-VortMag(i,j,k-1))/(2*theCellSize);
			//Gradient vector
			vec3 grad = vec3(gX, gY, gZ);

			//Normal (i,j,k)
			vec3 N = grad/(grad.Length() + 10.0e-20);

			//Omega
			vec3 omega = vec3(VortX(i,j,k), VortY(i,j,k), VortZ(i,j,k));

			//Force
			vec3 Fconf = ep*theCellSize*(N.Cross(omega));
			//Update Velocities with the vorticity confinement force
			target.mU(i,j,k) += dt*Fconf[0];
			target.mV(i,j,k) += dt*Fconf[1];
			target.mW(i,j,k) += dt*Fconf[2];
		}
	}	
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBuoyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	
	// TODO: Solve Ap = d for pressure.
	double deltaX = theCellSize;
	double deltaY = theCellSize;
	double deltaZ = theCellSize;
	double rho = 1.0;
	double vBoundary = 0.0;

	// 1. Construct d
	GridData d;
	d.initialize();
	//Multipler for the div*u
	double density = rho;//1.0;//mD(i, j, k);
	double multiplier = -(((deltaX*deltaX)*density)/dt);
	//Compute the Divergence Matrix that we need for solving for pressures
	FOR_EACH_CELL {
		vec3 pX = vec3(i,j,k);
		//Divergence of cell - fluid flowing out based on cell dims
		double divU, divV, divW;
		////Div-U
		if (i+1 == theDim[MACGrid::X]) {
			//Case for upper boundary & in cell
			divU = (vBoundary - mU(i,j,k)) / deltaX; 
		} else if ( i == 0) {
		//	//Case for lower boundary
			divU = (mU(i+1,j,k) - vBoundary) / deltaX;
		} else {
			//Case General X
			divU = (mU(i+1,j,k) - mU(i,j,k)) / deltaX; //old with first if
		}

		//Div-V
		if (j+1 == theDim[MACGrid::Y]) {
			//Case for Y Upper Boundary
			divV =  (vBoundary - mV(i,j,k)) / deltaY;
		} else if (j == 0 ) {
		//	//Case for Y Lower Boundary
			divV = (mV(i,j+1,k) - vBoundary) / deltaY;
		} else {
			//Case General Y
			divV = (mV(i,j+1,k) - mV(i,j,k)) / deltaY; 
		}

		//Div-W
		if (k+1 == theDim[MACGrid::Z]) {
			//Case for Z Upper Boundary
			divW = (vBoundary - mW(i,j,k)) / deltaZ;
		} else if (k == 0 ){
			//Case for Z Lower boundary
			divW = (mW(i,j,k+1) - vBoundary) / deltaZ;
		} else {
			//Case General Z
			divW = (mW(i,j,k+1) - mW(i,j,k)) / deltaZ; //old with first if 
		}

		//Update Divergence value
		double div = multiplier*(divU + divV + divW);
		d(i, j, k) = div;
	}


	// 2. Construct A
	GridDataMatrix A = this->AMatrix;

	// Solving vars
	int maxIterations = 400;
	double tolerance = 0.001;

	// 3. Solve for p - store in target.mP (pressure)
	conjugateGradient(A, target.mP, d, maxIterations, tolerance);
	//Init current pressure
	mP = target.mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	//Constant Multiplier
	double m = dt/rho;

	//Update vels using target.mP
	//u_n+1 = u_n -(dt*rho)*deltaP
	//cout<<currD<<endl;
	// X FACES
	FOR_EACH_FACE_X {  	
		//Make sure we aren't on boundaries
		if (i > 0 && i < theDim[MACGrid::X]) {
			//Update X-vel using Pressure gradient in X
			double currP = mP(i, j, k);
			double prevP = mP(i-1,j,k);

			double dPX = (m*(currP - prevP))/deltaX;
			target.mU(i, j, k) = mU(i,j,k) - dPX;
		}
		else {
			//Boundary CONDITIONS!
			target.mU(i,j,k) = 0;
		}
	}
	// Y FACES
	FOR_EACH_FACE_Y {
			
		if (j > 0 && j < theDim[MACGrid::Y]) {
			//Get Pressures for cell
			double currP = mP(i,j,k);
			double prevP = mP(i,j-1,k);
			
			//Udpdate Y-vel using Pressure Gradient in Y
			double dPY = (m*(currP - prevP))/deltaY; ;
			target.mV(i, j, k) = mV(i,j,k) - dPY;
		}
		else {
			//Boundary CONDITIONS!
			target.mV(i,j,k) = 0;
		}	

	}
	//Z FACES
	FOR_EACH_FACE_Z {
		if (k > 0 && k < theDim[MACGrid::Z]) {
			//Get Pressures for cell
			double currP = mP(i,j,k);//mP(i, j, k);
			double prevP = mP(i,j,k-1);
			
			//Update Z-vel using Pressure Gradient in Z
			double dPZ = (m*(currP - prevP))/deltaZ;
			target.mW(i, j, k) = mW(i,j,k) - dPZ;
		}
		else {
			//Boundary CONDITIONS!
			target.mW(i,j,k) = 0;
		}	
	}

	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;

	//Verify out divergence is zero
	//checkDivergence();
}

void MACGrid::checkDivergence() {
	//UPDATE IF YOU CHANGE DT!
	double dt = 0.04;
	double deltaX = theCellSize;
	double deltaY = theCellSize;
	double deltaZ = theCellSize;
	double rho = 1.0;
	double vBoundary = 0.0;

	// 1. Construct d
	GridData d;
	d.initialize();
	//Multipler for the div*u
	double density = rho;//1.0;//mD(i, j, k);
	double multiplier = -(((deltaX*deltaX)*density)/dt);
	//Compute the Divergence Matrix that we need for solving for pressures
	FOR_EACH_FACE {

		vec3 pX = vec3(i,j,k);
		//Divergence of cell - fluid flowing out based on cell dims
		double divU, divV, divW;
		//Div-U
		if (i == theDim[MACGrid::X]-1) {
			//Case for upper boundary & in cell
			divU = (vBoundary - mU(i,j,k)) / deltaX; 
		} else if ( i == 0) {
			//Case for lower boundary
			divU = (mU(i+1,j,k) - vBoundary) / deltaX;
		} else {
			//Case General X
			divU = (mU(i+1,j,k) - mU(i,j,k)) / deltaX; //old with first if
		}

		//Div-V
		if (j == theDim[MACGrid::Y]-1) {
			//Case for Y Upper Boundary
			divV =  (vBoundary - mV(i,j,k)) / deltaY;
		} else if (j == 0 ) {
			//Case for Y Lower Boundary
			divV = (mV(i,j+1,k) - vBoundary) / deltaY;
		}
		else {
			//Case General Y
			divV = (mV(i,j+1,k) - mV(i,j,k)) / deltaY; 
		}

		//Div-W
		if (k == theDim[MACGrid::Z]-1) {
			//Case for Z Upper Boundary
			divW = (vBoundary - mW(i,j,k)) / deltaZ;
		} else if (k == 0 ){
			//Case for Z Lower boundary
			divW = (mW(i,j,k+1) - vBoundary) / deltaZ;
		} else {
			//Case General Z
			divW = (mW(i,j,k+1) - mW(i,j,k)) / deltaZ; //old with first if 
		}

		//Update Divergence value
		double div = multiplier*(divU + divV + divW);
		d(i, j, k) = div;
	}

	FOR_EACH_FACE {
		//cout<<d(i,j,k)<<endl;
		if (d(i,j,k)> EPSILON) cout<< "INVALID DIVERGENCE!";
	}
}


vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}


/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(1.0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
