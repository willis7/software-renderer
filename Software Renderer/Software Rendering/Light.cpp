#include "Light.h"

Light::Light()
{
	ambient = Color(0.0f,0.0f,0.0f);
	diffuse = Color(1.0f,1.0f,1.0f);
	specular = Color(0.0f, 0.0f, 0.0f);
	worldVec = Vector3(0.0f, 0.0f, 0.0f);
	pointLight = false;
}

Light::Light(const Color &pAmbient, const Color &pDiffuse)
{
	ambient = pAmbient;
	diffuse = pDiffuse;
	specular = Color(0.0f, 0.0f, 0.0f);
	worldVec = Vector3(0.0f, 0.0f, 0.0f);
	pointLight = false;
}

void Light::setWorldVector(const Vector3 &vec)
{
	worldVec = vec;
}

bool Light::isPointLight()
{
	return pointLight;
}

void Light::setViewVector(const Vector3 &vec)
{
	viewVec = vec;
}

const Vector3 Light::worldVector()
{
	return worldVec;
}

const Vector3 Light::viewVector()
{
	return viewVec;
}

void Light::setAmbient(const Color &pAmbient)
{
	ambient = pAmbient;
}

void Light::setDiffuse(const Color &pDiffuse)
{
	diffuse = pDiffuse;
}	

void Light::setSpecular(const Color &pSpecular)
{
	specular = pSpecular;
}

const Color Light::getAmbient()
{
	return ambient;
}

const Color Light::getDiffuse()
{
	return diffuse;
}

const Color Light::getSpecular()
{
	return specular;
}
