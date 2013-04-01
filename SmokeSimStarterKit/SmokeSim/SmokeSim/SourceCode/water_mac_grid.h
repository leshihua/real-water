#ifndef WaterMACGrid_H_
#define WaterMACGrid_H_

#include "mac_grid.h"

class WaterMACGrid : public MACGrid {
public:
	// functions to override
	void reset();
	void advectTemperature(double dt);
	void advectDensity(double dt);

protected:
	// level set signed distance
	GridDataLSet mLSet;

};


#endif