#pragma once

#include "FrameBuffer.h"
#include "Matrix4.h"
#include "Light.h"

#include <vector>
using namespace std;


class RenderState
{
public:

	Light light;

	FrameBuffer fb;

	CullMode cullMode;			//None, Back Face or Front face

	RendererStyle style;		//Flat, Gouraud or Wireframe

	bool lightingEnabled;

	bool texturingEnabled;

	bool depthTest;

	int viewX;					//Viewport coordinates

	int viewY;

	Vector3 eyePos;				//Eye position

	Vector3 eyeDir;				//Eye direction

	//The stacks to hold both matrices
	vector<Matrix4> modelView;
	vector<Matrix4> projection;

	RenderState(void);
	
};
