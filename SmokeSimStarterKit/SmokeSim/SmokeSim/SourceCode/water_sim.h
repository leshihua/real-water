// ==========================================================================
// Copyright (C) 2009 Aline Normoyle
// Modified by Peter Kutz, 2011 and 2012.
// ==========================================================================
#ifndef WaterSim_H_
#define WaterSim_H_

#include "water_mac_grid.h"

class Camera;
class WaterSim
{
public:
   WaterSim();
   virtual ~WaterSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   virtual void setRecording(bool on, int width, int height);
   virtual bool isRecording();
	
	int getTotalFrames();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
	WaterMACGrid mGrid;
	bool mRecordEnabled;
	int mFrameNum;
	int mTotalFrameNum;
	
	int recordWidth;
	int recordHeight;
};

#endif