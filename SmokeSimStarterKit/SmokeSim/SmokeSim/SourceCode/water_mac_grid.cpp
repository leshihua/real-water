#include "water_mac_grid.h"

void WaterMACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);
   mLSet.initialize(0.0);
   setUpAMatrix();
}

void WaterMACGrid::advectTemperature(double dt) {
	assert("we don't advect temperature in a water sim!");
}

void WaterMACGrid::advectDensity(double dt) {
	assert("we don't advect density in a water sim!");
}
