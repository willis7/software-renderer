#include "FrameBuffer.h"

FrameBuffer::FrameBuffer() 
{
	width = 0;		height = 0;

	// allocate the buffer 
	colorBuffer = (u08 *)malloc(sizeof(u08) * width * height * 3);
	zBuffer = (float *)malloc(sizeof(float)* width * height);
}

FrameBuffer::~FrameBuffer(void)
{
	if (colorBuffer)
		free(colorBuffer);
	if (zBuffer)
		free(zBuffer);
}

void FrameBuffer::setSize(int w, int h)
{
	width = w;
	height = h;

	colorBuffer = (u08 *)malloc(sizeof(u08) * width * height * 3);
	zBuffer = (float *)malloc(sizeof(float)* width * height);

	clearColor();
	clearZ();
}

int FrameBuffer::getWidth(void)
{
	return width;
}

int FrameBuffer::getHeight(void)
{
	return height;
}

u08 *FrameBuffer::getColorPtr(int x, int y) 
{
	return (colorBuffer + ((y * width + x)*3));
}

void FrameBuffer::dumpToScreen(void) 
{
	glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, colorBuffer);
}

void FrameBuffer::clearColor()
{
	memset(colorBuffer, 0, sizeof(u08) * 3 * height * width);
}

void FrameBuffer::clearZ()
{
	memset(zBuffer, 0, sizeof(float) * height * width);
}

void FrameBuffer::putPixel(int x, int y,  float r, float g, float b)
{
	assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);

    r = (r < 0) ? 0 : 255 * ((r > 1) ? 1 : r);
    g = (g < 0) ? 0 : 255 * ((g > 1) ? 1 : g);
    b = (b < 0) ? 0 : 255 * ((b > 1) ? 1 : b);

    unsigned char *cp = &colorBuffer[3*(y * width + x)];

    cp[0] = (unsigned char) r;
    cp[1] = (unsigned char) g;
    cp[2] = (unsigned char) b;
}

void FrameBuffer::putPixel(int x, int y, float z, float r, float g, float b)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);

	int i = y*width + x;
	if (z >= zBuffer[i])
	{
	  r = (r < 0) ? 0 : 255 * ((r > 1) ? 1 : r);
	  g = (g < 0) ? 0 : 255 * ((g > 1) ? 1 : g);
	  b = (b < 0) ? 0 : 255 * ((b > 1) ? 1 : b);

	  unsigned char *cp = &colorBuffer[3*i];

	  cp[0] = (unsigned char) r;
	  cp[1] = (unsigned char) g;
	  cp[2] = (unsigned char) b;
	  zBuffer[i] = z;
	}
}
