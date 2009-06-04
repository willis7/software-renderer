#pragma once
#include "Globals.h"
#include "Vector3.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cmath>
#include <gl/glut.h>


class FrameBuffer
{
private:

	int width, height;

	/// all the buffers are defined from the lower left corner of the screen 
	u08 *colorBuffer;		// color buffer is unsigned bytes buffer size 3*w*h 
	float *zBuffer;

public:

	FrameBuffer();
	~FrameBuffer(void);

	void setSize(int w, int h);
	int getWidth(void);
	int getHeight(void);

	// get color ptr 
	u08 *getColorPtr(int x, int y);

	//These routines place a colour in the designated x,y coordinate
	void putPixel(int x, int y,  float r, float g, float b);
	void putPixel(int x, int y, float z,  float r, float g, float b);

	// this function dumps the data array to the screen using glDrawPixels()
	void dumpToScreen(void);

	//Clear the colorBuffer with all 0's
	void clearColor(void);
	void clearZ(void);
};
