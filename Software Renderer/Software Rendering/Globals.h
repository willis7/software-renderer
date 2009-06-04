#pragma once

typedef unsigned char u08;

#define PI 3.14159265359

//
struct Color
{
	float r, g, b;

	Color()
	{
		r = 0;
		g = 0;
		b = 0;
	}
	Color(float R, float G, float B)
	{
		r = R;
		g = G;
		b = B;
	}
};


enum CullMode
{
	CullNone,
	CullBack,
	CullFront
};

enum MatrixMode
{
	Projection,
	ModelView
};

enum ProjectionType
{
	Orthographic,
	Perspective
};

enum RendererStyle
{
	Flat = 0,
	Gouraud,
	Phong,
	Wireframe,
	NumRendererStyles
};
