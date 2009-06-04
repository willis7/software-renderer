#pragma once
#include "Vector3.h"
#include "Globals.h"

class Light
{
private:

	//Position in 3d space
	Vector3 worldVec;

	//view direction of light
	Vector3 viewVec;

	Color ambient;
	Color diffuse;
	Color specular;

	bool pointLight;

public:
	
	
	Light();
	Light(const Color &ambient, const Color &diffuse);

	bool isPointLight();
	void setWorldVector(const Vector3 &vec);	
	void setViewVector(const Vector3 &vec);
	const Vector3 worldVector();
	const Vector3 viewVector();
	void setAmbient(const Color &ambient);
	void setDiffuse(const Color &diffuse);
	void setSpecular(const Color &specular);
	const Color getAmbient();
	const Color getDiffuse();
	const Color getSpecular();
};
