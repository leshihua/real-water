#include "water_sim.h"
#include "constants.h"
#include "open_gl_headers.h"
#include "stb_image_write.h"
#include "custom_output.h"
#include "basic_math.h"

WaterSim::WaterSim() : mFrameNum(0), mTotalFrameNum(0), mRecordEnabled(false)
{
   reset();
}

WaterSim::~WaterSim() {
}

void WaterSim::reset()
{
   mGrid.reset();
	mTotalFrameNum = 0;
}

void WaterSim::step()
{

	double dt = 0.04; //0.1;
   // Step0: Gather user forces
   mGrid.updateSources();

   //TODO: Velocity extrapolation
   //TODO: Level Set Reinitialization

   // Step1: Calculate new velocities
   mGrid.advectSignedDistances(dt);
   mGrid.advectVelocity(dt);
   // Step2: Calculate new temperature
   mGrid.advectTemperature(dt);
   // Step3: Calculate new density 
   mGrid.advectDensity(dt);
   mGrid.addExternalForces(dt);
   //TODO: Level Set Remeshing
   mGrid.project(dt);
  
	
	mTotalFrameNum++;
}

void WaterSim::setRecording(bool on, int width, int height)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
	
	recordWidth = width;
	recordHeight = height;
}

bool WaterSim::isRecording()
{
   return mRecordEnabled;
}

void WaterSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) grabScreen();
}

void WaterSim::drawAxes()
{
	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
		glDisable(GL_LIGHTING);

		glLineWidth(2.0); 
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);
		glEnd();
	glPopAttrib();
}

void WaterSim::grabScreen()
{
	
	if (mFrameNum > 9999) exit(0);
	

	// TODO: Un-comment this to save each frame of smoke in CIS 460 volumetric text file format.
	/*
	// Save the smoke densities:
	char smoke_filename[2048];
	sprintf_s(smoke_filename, 2048, "smoke_%04d.txt", mFrameNum); // Use snprintf if your compiler supports C99.
	mGrid.saveSmoke(smoke_filename);
	*/
	

	// Save an image:

	unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

	for (int i=0; i<recordHeight; i++) 
	{
		glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
			bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
	}

	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "smoke_%04d.png", mFrameNum); 
	
	stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
	
	delete [] bitmapData;
	
	
	
	mFrameNum++;
	 
}

int WaterSim::getTotalFrames() {
	return mTotalFrameNum;
}
