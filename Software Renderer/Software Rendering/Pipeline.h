#pragma once
#include "RenderState.h"
#include "Light.h"


struct Texcoor{
	float u, v;
};

struct Triangle{
	Vector3 v[3]; 	//Original Verts
	Vector3 vt[3];	//Translated Verts
	Vector3 n;		//Normal
	Color c[3];	
	Texcoor	t[3];

	//Material properties
	Color specularcolor;
    Color diffusecolor;
    Color ambientcolor;
	float shininess;
};

class Pipeline 
{

private:

	RenderState state;

	static const float epsilon;

public:

	//Default constructor
	Pipeline();

//Buffer Operations
	void  draw();
	void  clear();

//Perspective
	void  lookAt(const Vector3 &eye, const Vector3 &center, const Vector3 &up);
	void  perspective(float fov, float aspect, float zNear, float zFar);

	void  setViewport(int x, int y, int width, int height);
	void  viewTransform(Vector3 &v);

	void perspectiveDivide(Vector3 &v, int w);

//Matrix operations
	void  popMatrix(MatrixMode mode);
	void  pushMatrix(MatrixMode mode);

	void  loadIdentity(MatrixMode mode);
	const Matrix4 getModelView();
	const Matrix4 getProjectionView();

	void  translate(const Vector3 &t);
	void  scale(float s);
	void  rotate(const Vector3 &axis, int angle);

//Depth testing and culling
	CullMode getCullMode();
	void  setCullMode(CullMode cullMode);
	void  enableDepthTest(bool enable);
	bool  depthTestEnabled();

	bool cullTriangle(Triangle &t);
	
//Texturing
	bool texturingEnabled();
	void enableTexturing(bool texturing);

//Lighting
	bool lightingEnabled();
	void enableLighting(bool lighting);

	const Vector3 getLightVector();
	void  setLightVector(const Vector3 &lightvec);

	const Color getLightAmbient();
	void  setLightAmbient(const Color &ambient);

	const Color getLightDiffuse();
	void  setLightDiffuse(const Color &diffuse);

	const Color getLightSpecular();
	void  setLightSpecular(const Color &specular);

	void calcLight(Triangle &t, vector<Light> &lights);


//Set the render style (Wireframe, Solid, Gouraud)

	void setStyle(RendererStyle style);

//Create Scene

	void sortVerts(Triangle &t);
	void putPixel(Vector3 p, Color c);
	void putPixel(float x, float y, float r, float g, float b);
	void drawLine(float r1, float g1, float b1, float x1, float y1, float r2, float g2, float b2, float x2, float y2);
	bool clipLine(float &x1, float &y1, float &x2, float &y2);
	void drawTriangle(Triangle &t, vector<Light> &lights);
	void drawTriangle(Triangle &t);

	void bezierCurve(Vector3 *points, int order, int LOD);

	void dump();

};

