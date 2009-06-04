#include <iostream>
#include "Pipeline.h"
#include <gl/glut.h>

using namespace std;

#define WINDOW_WIDTH 1000
#define WINDOW_HEIGHT 800

Pipeline pl;

//Keep an array of the lights in the scene
vector<Light> sceneLights;
Light light1;

//Points for bezier curve testing
Vector3 points[4] = {Vector3(700,700,0),
					 Vector3(300,400,100),
					 Vector3(500, 0,0),
					 Vector3(100,100,40)};


// Initializes the 3D rendering
void initRendering()
{	
	//Clear both buffers
	pl.clear();

	//Enables the lighting calculations
	pl.enableLighting(true);

	light1.setDiffuse(Color(1.0, 1.0, 1.0));
	light1.setAmbient(Color(1.0, 0.0, 0.0));
	light1.setSpecular(Color(1.0, 1.0, 1.0));
	light1.setWorldVector(Vector3(10.0,0.0,-1.0));
	
	sceneLights.push_back(light1);
	
	//To test hidden surface removal the parameter can be changed to CullFront
	pl.setCullMode(CullBack);																	
}

//Called when window is resized
void handleResize(int w, int h)			
{
	pl.clear();
	pl.setViewport(0, 0, w, h);

	pl.loadIdentity(Projection);

	pl.perspective(60.0,				//Camera angle	
					w/h,				// The Width-Height ratio
					-1.0,				// The near Z clipping co-ordinate (dictates how close to the eye it is allowed to come before not getting shown)
					-1000.0);			// The far Z clipping co-ordinates (dictates how far from the eye it is allowed to go before no longer getting shown)
}

/////////////////////////////////////////////////////////
//////////// MAIN DRAWING ROUTINE ///////////////////////
/////////////////////////////////////////////////////////
void drawScene()
{


// TRIANGLE 1
	Triangle t;

	//Vertex1
	t.c[0].r = 1.0; 	
	t.c[0].g = 0.0;	
	t.c[0].b = 0.0;	
	//Vertex2
	t.c[1].r = 0.0;
	t.c[1].g = 1.0; 
	t.c[1].b = 0.0;
	//Vertex3
	t.c[2].r = 0.0;
	t.c[2].g = 0.0;
	t.c[2].b = 1.0;

	t.n = Vector3(0.0, 0.0, -1.0);

	t.ambientcolor.r = 1.0; t.diffusecolor.r = 0.5; t.specularcolor.r = 0.992157;
	t.ambientcolor.g = 0.0; t.diffusecolor.g = 0.5; t.specularcolor.g = 0.941176;
	t.ambientcolor.b = 0.0; t.diffusecolor.b = 0.5; t.specularcolor.b = 0.807843;

	t.shininess = 1;

	t.v[0].x = 100; t.v[1].x = 200;	t.v[2].x = 100;
	t.v[0].y = 200;	t.v[1].y = 200; t.v[2].y = 100;
	t.v[0].z = 10;	t.v[1].z = 10;	t.v[2].z = 10;	
	

	pl.setStyle(Gouraud);			
	pl.drawTriangle(t, sceneLights);



//  TRIANGLE 2
	Triangle t1;
	//Vertex1
	t1.c[0].r = 0.0; 	
	t1.c[0].g = 1.0;
	t1.c[0].b = 0.0; 
	//Vertex2
	t1.c[1].r = 1.0;  
	t1.c[1].g = 0.0;
	t1.c[1].b = 0.0;
	//Vertex3
	t1.c[2].r = 0;
	t1.c[2].g = 0;	
	t1.c[2].b = 1.0;

	t1.n = Vector3(0.0, 0.0, -1.0);

	t1.ambientcolor.r = 1.0; t1.diffusecolor.r = 0.5;
	t1.ambientcolor.g = 0.0; t1.diffusecolor.g = 0.5;
	t1.ambientcolor.b = 0.0; t1.diffusecolor.b = 0.5;

	t1.shininess = 1;

	t1.v[0].x = 150;
	t1.v[0].y = 200;
	t1.v[0].z = 100; 

	t1.v[1].x = 150;
	t1.v[1].y = 100;
	t1.v[1].z = 100; 

	t1.v[2].x = 50;	
	t1.v[2].y = 100;	 
	t1.v[2].z = 100;	
	

	pl.setStyle(Gouraud);	
	pl.translate(Vector3(50,0,0));						//Translates the second triangle in the x by 50
	pl.drawTriangle(t1, sceneLights);





// This calls for opengl to dump my frame buffer and then refresh the screen
	pl.dump();
	glutSwapBuffers();
	glFlush();
}

int main(int argc, char** argv)
{

	//Initialise glut
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);						//Set the window size

	glutCreateWindow("Software Renderer");
	initRendering();														//Initialise rendering

	//Set handler functions for drawing, keypresses and window resizes
	glutDisplayFunc(drawScene);
	glutReshapeFunc(handleResize);

	glutMainLoop();															//Start the main loop..... glutMainLoop never returns!! 
																			//(tells glut to take control of the program i.e. handle the resizes)

	return 0;																// This is never reached!!
}