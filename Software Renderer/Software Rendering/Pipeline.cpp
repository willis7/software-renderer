#include "Pipeline.h"

const float Pipeline::epsilon = 0.0001f;



Pipeline::Pipeline()
{
	//Set the top matrix of both modelview and projection to the Id matrix
	state.modelView.push_back(Matrix4());
	state.projection.push_back(Matrix4());	
}

void Pipeline::loadIdentity(MatrixMode mode)
{
	if (mode == ModelView)
    {
		state.modelView.back().identity();  
    }
    else
    {
      state.projection.back().identity();
    }
}

void Pipeline::popMatrix(MatrixMode mode)
{
	if(mode == ModelView && state.modelView.size() > 1)
	{
		state.modelView.pop_back();

	}else if(mode == Projection && state.projection.size() > 1)
	{
		state.projection.pop_back();
	}
}

void Pipeline::pushMatrix(MatrixMode mode)
{
	if (mode == ModelView)
		state.modelView.push_back(state.modelView.back());
	else if (mode == Projection)
		state.projection.push_back(state.projection.back());
 
}

const Matrix4 Pipeline::getModelView()
{
	return state.modelView.back();
}

const Matrix4 Pipeline::getProjectionView()
{
	return state.projection.back();
}

void Pipeline::clear()
{
	state.fb.clearColor();
	state.fb.clearZ();
}

void Pipeline::lookAt(const Vector3 &eye, const Vector3 &center, const Vector3 &up)
{
	//eye 		Specifies the position of the eye point.
	//center 	Specifies the position of the reference point.
	//up 		Specifies the direction of the up vector.

	Vector3 f(center - eye);
	f.normalize();

	Vector3 u(up);
	u.normalize();

	Vector3 s(crossProd(f, u));
	Vector3 v(crossProd(s, f));

	Matrix4 m;
		
	m.mat4[0] = s.x;  m.mat4[4] = s.y;  m.mat4[8] = s.z;
	m.mat4[1] = v.x;  m.mat4[5] = v.y;  m.mat4[9] = v.z;
	m.mat4[3] = -f.x; m.mat4[6] = -f.y; m.mat4[10] = -f.z;

	//Set this to the current matrix (top of the stack)
	state.modelView.back() = state.modelView.back() * m;
	
}

void Pipeline::perspective(float fov, float aspect, float zNear, float zFar)
{

	Matrix4 m;
	float f = 1/tan(fov);

	m.mat4[0] = f/aspect;
	m.mat4[5] = f;
	m.mat4[10] = (zFar + zNear)/(zFar - zNear);
	m.mat4[11] = -1;
	m.mat4[14] = (2*zFar*zNear )/(zNear - zFar);
	m.mat4[15] = 0;	

	state.projection.back() = state.projection.back() * m;	

}

void  Pipeline::setViewport(int x, int y, int width, int height)
{
	state.viewX = x;
	state.viewY = y;
	state.fb.setSize(width, height);
}

void  Pipeline::viewTransform(Vector3 &v)
{
	v.x = (v.x + 1)*(state.fb.getWidth() / 2) + state.viewX;
	v.y = (v.y + 1)*(state.fb.getHeight() / 2) + state.viewY;
}

void Pipeline::perspectiveDivide(Vector3 &v, int w)
{
	v.x = v.x/w;
	v.y = v.y/w;
	v.z = v.z/w;

}

void  Pipeline::translate(const Vector3 &t)
{
	Matrix4 m, mv;

	m.mat4[12] = t.x;
	m.mat4[13] = t.y;
	m.mat4[14] = t.z;

	//Set this to the current matrix (top of the stack)
	state.modelView.back() = state.modelView.back() * m;
	
}

void  Pipeline::scale(float s)
{
	Matrix4 m;

	m.mat4[0] = s;
	m.mat4[5] = s;
	m.mat4[10] = s;

	//Set this to the current matrix (top of the stack)
	state.modelView.back() = state.modelView.back() * m;
}

void  Pipeline::rotate(const Vector3 &axis, int a)
{
	Matrix4 m;

	if(axis.x == 1)
	{
		m.mat4[5] = Matrix4::COS_LUT[a];
		m.mat4[6] = Matrix4::SIN_LUT[a];
		m.mat4[9] = -Matrix4::SIN_LUT[a];
		m.mat4[10] = Matrix4::COS_LUT[a];

	}else if(axis.y == 1)
	{
		m.mat4[0] = Matrix4::COS_LUT[a];
		m.mat4[2] = -Matrix4::SIN_LUT[a];
		m.mat4[8] = Matrix4::SIN_LUT[a];
		m.mat4[10] = Matrix4::COS_LUT[a];

	}else if(axis.z == 1)
	{
		m.mat4[0] = Matrix4::COS_LUT[a];
		m.mat4[1] = Matrix4::SIN_LUT[a];
		m.mat4[4] = -Matrix4::SIN_LUT[a];
		m.mat4[5] = Matrix4::COS_LUT[a];
	}else
	{
		printf("ERROR: no axis was selected\n");
	}

	//Set this to the current matrix (top of the stack)
	state.modelView.back() = state.modelView.back() * m;

}
CullMode Pipeline::getCullMode()
{
	return state.cullMode;
}

void  Pipeline::setCullMode(CullMode cullMode)
{
	state.cullMode = cullMode;
}

void  Pipeline::enableDepthTest(bool enable)
{
	state.depthTest = enable;
}

bool  Pipeline::depthTestEnabled()
{
	return state.depthTest;
}

bool Pipeline::texturingEnabled()
{
	return state.texturingEnabled;
}

void Pipeline::enableTexturing(bool texturing)
{
	state.texturingEnabled = texturing;
}

bool Pipeline::lightingEnabled()
{
	return state.lightingEnabled;
}

void Pipeline::enableLighting(bool lighting)
{
	state.lightingEnabled = lighting;
}

const Vector3 Pipeline::getLightVector()
{
	return state.light.worldVector();
}

void  Pipeline::setLightVector(const Vector3 &lightvec)
{
	state.light.setWorldVector(lightvec);

	Vector3 v =  state.modelView.back()* lightvec;

	if(lightvec.z < epsilon)
		v.normalize();
	state.light.setViewVector(v);
}

const Color Pipeline::getLightAmbient()
{
	return state.light.getAmbient();
}

void  Pipeline::setLightAmbient(const Color &ambient)
{
	state.light.setAmbient(ambient);
}

const Color Pipeline::getLightDiffuse()
{
	return state.light.getDiffuse();
}

void  Pipeline::setLightDiffuse(const Color &diffuse)
{
	return state.light.setDiffuse(diffuse);
}

const Color Pipeline::getLightSpecular()
{
	return state.light.getSpecular();
}

void  Pipeline::setLightSpecular(const Color &specular)
{
	state.light.setSpecular(specular);
}
void Pipeline::calcLight(Triangle &t, vector<Light> &lights)
{
	//Color = Ka*ambientColor + SumOf(Kd*diffuseColor*(N dot L) + Ks*specularColor*(R dot V)^shininess)

	Color ka = t.ambientcolor;
	Color kd = t.diffusecolor;
	Color ks = t.specularcolor;

	//ambientColor = white
	Color ambientColor(1.0, 1.0, 1.0);

	//Ka*ambientColor
	float ambientFactorR = ka.r*ambientColor.r;
	float ambientFactorG = ka.g*ambientColor.g;
	float ambientFactorB = ka.b*ambientColor.b;

	//Kd*diffuseColor*(N dot L)
	float diffuseFactorR, diffuseFactorG, diffuseFactorB;
	diffuseFactorR = diffuseFactorG = diffuseFactorB = 0;

	//Ks*specularColor*(R dot V)^shininess
	float specularFactorR, specularFactorG, specularFactorB;
	specularFactorR = specularFactorG = specularFactorB = 0;

	//for each pixel in the triangle
	for(int i = 0; i < 3; i++)
	{
		//for all the lights in the scene
		for(int j = 0; j < lights.size(); j++)
		{
			
			//diffuseColor = surfaceColor * lightColor
			//Color diffuseColor(lights[j].getDiffuse().r *kd.r, lights[j].getDiffuse().g*kd.g, lights[j].getDiffuse().b*kd.b);

			//specularColor = lightColor
			Color specularColor(lights[j].getSpecular());

			//Kd*diffuseColor*(N dot L)
			float dotProd;
			dotProd = t.n * lights[j].worldVector();

			diffuseFactorR += kd.r * lights[j].getDiffuse().r * dotProd;
			diffuseFactorG += kd.g * lights[j].getDiffuse().g * dotProd;
			diffuseFactorB += kd.b * lights[j].getDiffuse().b * dotProd;

			//Ks*specularColor*(R dot V)^shininess
			//R = 2*N*(N dot L) - L
			Vector3 R;
			R = 2*t.n*dotProd - lights[j].worldVector();
			float dotProd2;
			dotProd2 = R * state.eyePos;			

			specularFactorR += ka.r * specularColor.r * dotProd2;
			specularFactorG += ka.g * specularColor.g * dotProd2;
			specularFactorB += ka.b * specularColor.b * dotProd2;

		}

		//Finally set the triangles color values
		t.c[i].r = ambientFactorR + diffuseFactorR + pow(specularFactorR, t.shininess);
		t.c[i].g = ambientFactorG + diffuseFactorG + pow(specularFactorG, t.shininess);
		t.c[i].b = ambientFactorB + diffuseFactorB + pow(specularFactorB, t.shininess);

		//Clear ready for next vert
		diffuseFactorR = diffuseFactorG = diffuseFactorB = 0;
		specularFactorR = specularFactorG = specularFactorB = 0;
	}


}
void Pipeline::setStyle(RendererStyle s)
{
	state.style = s;
}

void Pipeline::dump()
{
	state.fb.dumpToScreen();
}

void Pipeline::sortVerts(Triangle &t)
{
	Vector3 temp1;

	//v[0].y <= v[1].y <= v[2].y
	if( t.vt[2].y < t.vt[1].y)
	{
		temp1 = t.vt[1];
		t.vt[1] = t.vt[2];
		t.vt[2] = temp1;
	}
	if( t.vt[1].y < t.vt[0].y)
	{
		temp1 = t.vt[0];
		t.vt[0] = t.vt[1];
		t.vt[1] = temp1;
	}
	if(t.vt[2].y < t.vt[1].y)
	{
		temp1 = t.vt[1];
		t.vt[1] = t.vt[2];
		t.vt[2] = temp1;
	}
}

void Pipeline::putPixel(Vector3 p, Color c)
{
	state.fb.putPixel(p.x, p.y, p.z, c.r, c.g, c.b);
}

void Pipeline::putPixel(float x, float y, float r, float g, float b)
{
	state.fb.putPixel(x, y, r, g, b);
}

void Pipeline::drawLine(float r1, float g1, float b1, float startx, float starty, float r2, float g2, float b2, float endx, float endy)
{
	// We don't want to alter the original data
	float x1 = startx;
	float y1 = starty;
	float x2 = endx;
	float y2 = endy;

	// Clip the line to fall within the screen
	if(clipLine(x1, y1, x2, y1) == false)
		return;


	float xdiff = (x2 - x1);
    float ydiff = (y2 - y1);

    if(xdiff == 0.0f && ydiff == 0.0f) 
	{
		putPixel(x1, y1, r1, g1, b1);
        return;
    }

	if(fabs(xdiff) > fabs(ydiff)) 
	{
		float xmin, xmax;

		// set xmin to the lower x value given
		// and xmax to the higher value
		if(x1 < x2) 
		{
			xmin = x1;
			xmax = x2;
		} else {
			xmin = x2;
			xmax = x1;
		}

		// draw line in terms of y slope
		float slope = ydiff / xdiff;
		for(float x = xmin; x <= xmax; x += 1.0f) 
		{

			float y = y1 + ((x - x1) * slope);
            
			//Interpolate the color
			float rNew	= r1 + ((r2 - r1) * ((x - x1) / xdiff));
			float gNew	= g1 + ((g2 - g1) * ((x - x1) / xdiff));
			float bNew	= b1 + ((b2 - b1) * ((x - x1) / xdiff));

			putPixel(x, y, rNew, gNew, bNew);
		}

	}else 
	{

        float ymin, ymax;

        // set ymin to the lower y value given
        // and ymax to the higher value
        if(y1 < y2) 
		{
                ymin = y1;
                ymax = y2;
        } else 
		{
                ymin = y2;
                ymax = y1;
        }

        // draw line in terms of x slope
        float slope = xdiff / ydiff;
        for(float y = ymin; y <= ymax; y += 1.0f) 
		{
                float x = x1 + ((y - y1) * slope);
                
				//Interpolate the color
				float rNew	= r1 + ((r2 - r1) * ((y - y1) / ydiff));
				float gNew	= g1 + ((g2 - g1) * ((y - y1) / ydiff));
				float bNew	= b1 + ((b2 - b1) * ((y - y1) / ydiff));

                putPixel(x, y, rNew, gNew, bNew);
        }
	}
}

bool Pipeline::clipLine(float &x1, float &y1, float &x2, float &y2)
{
	double slope =	(double)(y2 - y1) /
					(double)(x2 - x1);

	int clipS = 0, clipE = 0;

	// Set resolution = frame buffer
	int resWidth = state.fb.getWidth();
	int resHeight = state.fb.getHeight();


	// Loop while at least one point is outside the canvas.
	do
	{
		// Location tests for the start point.
		clipS = (( x1 < 0) << 3) | ((x1 >= resWidth) << 2) |
		((y1 < 0) << 1) | (y1 >= resHeight);

		// Location tests for the end point.
		clipE = ((x2 < 0) << 3) | ((x2 >= resWidth) << 2) |
		((y2 < 0) << 1) | (y2 >= resHeight);

		// The line is completely outside of the canvas.
		if(clipS & clipE)
			return false;

		// If the start is outside then clip it based on location.
		if(clipS)
		{
			// Left side.
			if(clipS & 8)
			{
			y1 -= (int)((double)x1 * slope);
			x1 = 0;
			}
			else
			{
				// Right side.
				if(clipS & 4)
				{
				y1 += (int)((double)(resWidth - x1) *
				slope);
				x1 = resWidth - 1;
				}
				else
				{
					// Top side.
					if(clipS & 2)
					{
					x1 -= (int)((double)y1 / slope);
					y1 = 0;
				}
				else
				{
					// Bottom side.
					if(clipS & 1)
					{
						x1 += (int)((double)(resHeight - y1) / slope);

						y1 = resHeight - 1;
					}
				}
			}
		}
	}

	// Clip end point if it is anywhere outside of the canvas.
	if(clipE)
	{
		// Left side.
		if(clipE & 8)
		{
			y2 += (int)((double)(0 - x2) * slope);
			x2 = 0;
		}
		else
		{
			// Right side.
			if(clipE & 4)
			{
				y2 += (int)((double)(resWidth - x2) * slope);

				x2 = resWidth - 1;
			}
			else
			{
				// Top side.
				if(clipE & 2)
				{
					x2 += (int)((double)(0 - y2) / slope);
					y2 = 0;
				}
				else
				{
					// Bottom side.
					if(clipE & 1)
					{
						x2 += (int)((double)(resHeight - y2) / slope);
						y2 = resHeight - 1;
					}
				}
			}
		}
	}
	} while(clipS | clipE);
	return true;

}
bool Pipeline::cullTriangle(Triangle &t)
{
	// Check if the triangle already has a normal
	if(t.n == Vector3(0.0, 0.0, 0.0))
	{
		//compute face normal from the points
		// counter clockwise ordering is assumed
		Vector3 v = t.v[0] - t.v[1];
		Vector3 w = t.v[2] - t.v[1];

		t.n = crossProd(v, w);
		t.n.normalize();
	}


	if(state.cullMode == 1)			//Back face cull
	{
		//compute the see which way the face is facing in comparison
		//  to the camera (located at 0,0,0)
		float angle = state.eyeDir * t.n;

		if( angle < 0)
			return true;
		else
			return false;
	}
	else							//Front face cull
	{
		//compute the see which way the face is facing in comparison
		//  to the camera (located at 0,0,0)
		float angle = state.eyeDir * t.n;

		if( angle > 0)
			return true;
		else
			return false;
	}
}
void Pipeline::drawTriangle(Triangle &t, vector<Light> &lights)
{
	Matrix4 mv;

	mv = getModelView();

	//Clear previous values
	t.vt[0].zero();
	t.vt[1].zero();
	t.vt[2].zero();

	// MULTIPLY MODELVIEW
	t.vt[0] = mv * t.v[0];
	t.vt[1] = mv * t.v[1];
	t.vt[2] = mv * t.v[2];

	mv = getProjectionView();
	
	// MULTIPLY PROJECTION
	t.vt[0] = mv * t.vt[0];
	t.vt[1] = mv * t.vt[1];
	t.vt[2] = mv * t.vt[2];

	// PERSPECTIVE DIVIDE
	perspectiveDivide(t.vt[0], 1);
	perspectiveDivide(t.vt[1], 1);
	perspectiveDivide(t.vt[2], 1);

	// VIEWPORT TRANSFORM
	//viewTransform(t.vt[0]);
	//viewTransform(t.vt[1]);
	//viewTransform(t.vt[2]);

	//Test for horizontal and verticle lines
	if((t.vt[0].x == t.vt[1].x && t.vt[1].x == t.vt[2].x) || (t.vt[0].y == t.vt[1].y && t.vt[1].y == t.vt[2].y))
		return;

	//CULLING
	//Check which way the triangle is facing and depending on mode
	//cull accordingly
	if(state.cullMode != CullNone)
	{
		if(cullTriangle(t) == false)
			return;
	}

	//SORT TRIANGLE VERTS FOR TRIANGLE DRAW ALGORITHM
	sortVerts(t);

	//LIGHTING
	if(state.lightingEnabled == true)
	{
		calcLight(t, lights);
	}

//DRAW FLAT
	if(state.style == Flat)
	{

		float dx1, dx2, dx3;
		
		if (t.vt[1].y - t.vt[0].y > 0) 
			dx1 = (t.vt[1].x - t.vt[0].x)/(t.vt[1].y - t.vt[0].y);
		else 
			dx1 = t.vt[1].x - t.vt[0].x;

		if (t.vt[2].y - t.vt[0].y > 0) 
			dx2 = (t.vt[2].x - t.vt[0].x)/(t.vt[2].y - t.vt[0].y);
		else 
			dx2 = 0;

		if (t.vt[2].y - t.vt[1].y > 0) 
			dx3 = (t.vt[2].x - t.vt[1].x)/(t.vt[2].y - t.vt[1].y);
		else 
			dx3 = 0;

		Vector3 S = t.vt[0];
		Vector3 xe = t.vt[0];

		if(dx1 > dx2) 
		{
			for( ; S.y <= t.vt[1].y; S.y++, xe.y++, S.x += dx2, xe.x += dx1)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			xe = t.vt[1];

			for( ; S.y <= t.vt[2].y; S.y++, xe.y++, S.x += dx2, xe.x += dx3)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			
		}else
		{
			for( ; S.y <= t.vt[1].y; S.y++, xe.y++, S.x += dx1, xe.x += dx2)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			S = t.vt[1];

			for( ; S.y <= t.vt[2].y; S.y++, xe.y++, S.x += dx3, xe.x += dx2)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}
		}				
	}	

//DRAW GOURAUD
	if(state.style == Gouraud)
	{
		
		//deltas used in interpolation of x-coordinate 
		float dx1,dx2,dx3;
		//deltas used in interpolation of color's components
		float dr1,dr2,dr3, dg1,dg2,dg3, db1,db2,db3;

		//Initialise the Delta's

		// (B.y - A.y > 0)
		if (t.vt[1].y - t.vt[0].y > 0) 
		{
			dx1 = (t.vt[1].x - t.vt[0].x)/(t.vt[1].y - t.vt[0].y);
			dr1 = (t.c[1].r - t.c[0].r)/(t.vt[1].y - t.vt[0].y); 
			dg1 = (t.c[1].g - t.c[0].g)/(t.vt[1].y - t.vt[0].y); 
			db1 = (t.c[1].b - t.c[0].b)/(t.vt[1].y - t.vt[0].y); 
		} else 
		{
			dx1 = t.vt[1].x - t.vt[0].x;
			dr1 = t.c[1].r - t.c[0].r;
			dg1 = t.c[1].g - t.c[0].g;
			db1 = t.c[1].b - t.c[0].b;
		}

		// (C.y - A.y > 0)
		if (t.vt[2].y - t.vt[0].y > 0) 
		{
			dx2=(t.vt[2].x - t.vt[0].x)/(t.vt[2].y - t.vt[0].y);
			dr2=(t.c[2].r - t.c[0].r)/(t.vt[2].y - t.vt[0].y);
			dg2=(t.c[2].g - t.c[0].g)/(t.vt[2].y - t.vt[0].y);
			db2=(t.c[2].b - t.c[0].b)/(t.vt[2].y - t.vt[0].y);
		} else 
			dx2=dr2=dg2=db2=0;

		// (C.y - B.y > 0)
		if (t.vt[2].y - t.vt[1].y > 0) 
		{
			dx3=(t.vt[2].x - t.vt[1].x)/(t.vt[2].y - t.vt[1].y);
			dr3=(t.c[2].r - t.c[1].r)/(t.vt[2].y - t.vt[1].y);
			dg3=(t.c[2].g - t.c[1].g)/(t.vt[2].y - t.vt[1].y);
			db3=(t.c[2].b - t.c[1].b)/(t.vt[2].y - t.vt[1].y);
		} else 
			dx3=dr3=dg3=db3=0;

		Vector3 S = t.vt[0];
		Vector3 E = t.vt[0];

		Color cs = t.c[0];
		Color ce = t.c[0];

		float dr, dg, db;
		Vector3 p;
		Color c;

		if(dx1 > dx2) 
		{
			for(; S.y <= t.vt[1].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;
				
				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx2; cs.r += dr2; cs.g += dg2; cs.b += db2;
				E.x += dx1; ce.r += dr1; ce.g += dg1; ce.b += db1;
			}

			E = t.vt[1]; 
			ce = t.c[1];

			for(; S.y <= t.vt[2].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx2; cs.r += dr2; cs.g += dg2; cs.b += db2;
				E.x += dx3; ce.r += dr3; ce.g += dg3; ce.b += db3;
			}

		}else
		{
			for(; S.y <= t.vt[1].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx1; cs.r += dr1; cs.g += dg1; cs.b += db1;
				E.x += dx2; ce.r += dr2; ce.g += dg2; ce.b += db2;
			}

			S = t.vt[1]; 
			cs = t.c[1];

			for(; S.y <= t.vt[2].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx3; cs.r += dr3; cs.g += dg3; cs.b += db3;
				E.x += dx2; ce.r += dr2; ce.g += dg2; ce.b += db2;
			}
		}
	}

//DRAW WIREFRAME
	if(state.style == Wireframe)
	{
		drawLine( t.c[0].r, t.c[0].g, t.c[0].b, t.vt[0].x, t.vt[0].y, t.c[1].r, t.c[1].g, t.c[1].b, t.vt[1].x, t.vt[1].y );
		drawLine( t.c[1].r, t.c[1].g, t.c[1].b, t.vt[1].x, t.vt[1].y, t.c[2].r, t.c[2].g, t.c[2].b, t.vt[2].x, t.vt[2].y );
		drawLine( t.c[0].r, t.c[0].g, t.c[0].b, t.vt[0].x, t.vt[0].y, t.c[2].r, t.c[2].g, t.c[2].b, t.vt[2].x, t.vt[2].y );		
	}

}

void Pipeline::drawTriangle(Triangle &t)
{
	Matrix4 mv;

	mv = getModelView();

	//Clear previous values
	t.vt[0].zero();
	t.vt[1].zero();
	t.vt[2].zero();

	// MULTIPLY MODELVIEW
	t.vt[0] = mv * t.v[0];
	t.vt[1] = mv * t.v[1];
	t.vt[2] = mv * t.v[2];

	mv = getProjectionView();
	
	// MULTIPLY PROJECTION
	t.vt[0] = mv * t.vt[0];
	t.vt[1] = mv * t.vt[1];
	t.vt[2] = mv * t.vt[2];

	// PERSPECTIVE DIVIDE
	perspectiveDivide(t.vt[0], 1);
	perspectiveDivide(t.vt[1], 1);
	perspectiveDivide(t.vt[2], 1);

	// VIEWPORT TRANSFORM
	//viewTransform(t.vt[0]);
	//viewTransform(t.vt[1]);
	//viewTransform(t.vt[2]);

	//Test for horizontal and verticle lines
	if((t.vt[0].x == t.vt[1].x && t.vt[1].x == t.vt[2].x) || (t.vt[0].y == t.vt[1].y && t.vt[1].y == t.vt[2].y))
		return;

	//CULLING
	//Check which way the triangle is facing and depending on mode
	//cull accordingly
	if(state.cullMode != CullNone)
	{
		if(cullTriangle(t) == false)
			return;
	}

	//SORT TRIANGLE VERTS FOR TRIANGLE DRAW ALGORITHM
	sortVerts(t);


//DRAW FLAT
	if(state.style == Flat)
	{

		float dx1, dx2, dx3;
		
		if (t.vt[1].y - t.vt[0].y > 0) 
			dx1 = (t.vt[1].x - t.vt[0].x)/(t.vt[1].y - t.vt[0].y);
		else 
			dx1 = t.vt[1].x - t.vt[0].x;

		if (t.vt[2].y - t.vt[0].y > 0) 
			dx2 = (t.vt[2].x - t.vt[0].x)/(t.vt[2].y - t.vt[0].y);
		else 
			dx2 = 0;

		if (t.vt[2].y - t.vt[1].y > 0) 
			dx3 = (t.vt[2].x - t.vt[1].x)/(t.vt[2].y - t.vt[1].y);
		else 
			dx3 = 0;

		Vector3 S = t.vt[0];
		Vector3 xe = t.vt[0];

		if(dx1 > dx2) 
		{
			for( ; S.y <= t.vt[1].y; S.y++, xe.y++, S.x += dx2, xe.x += dx1)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			xe = t.vt[1];

			for( ; S.y <= t.vt[2].y; S.y++, xe.y++, S.x += dx2, xe.x += dx3)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			
		}else
		{
			for( ; S.y <= t.vt[1].y; S.y++, xe.y++, S.x += dx1, xe.x += dx2)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}

			S = t.vt[1];

			for( ; S.y <= t.vt[2].y; S.y++, xe.y++, S.x += dx3, xe.x += dx2)
			{
				drawLine(t.c[0].r, t.c[0].g, t.c[0].b, S.x, S.y,
							t.c[0].r, t.c[0].g, t.c[0].b, xe.x, S.y);
			}
		}				
	}	

//DRAW GOURAUD
	if(state.style == Gouraud)
	{
		
		//deltas used in interpolation of x-coordinate 
		float dx1,dx2,dx3;
		//deltas used in interpolation of color's components
		float dr1,dr2,dr3, dg1,dg2,dg3, db1,db2,db3;

		//Initialise the Delta's

		// (B.y - A.y > 0)
		if (t.vt[1].y - t.vt[0].y > 0) 
		{
			dx1 = (t.vt[1].x - t.vt[0].x)/(t.vt[1].y - t.vt[0].y);
			dr1 = (t.c[1].r - t.c[0].r)/(t.vt[1].y - t.vt[0].y); 
			dg1 = (t.c[1].g - t.c[0].g)/(t.vt[1].y - t.vt[0].y); 
			db1 = (t.c[1].b - t.c[0].b)/(t.vt[1].y - t.vt[0].y); 
		} else 
		{
			dx1 = t.vt[1].x - t.vt[0].x;
			dr1 = t.c[1].r - t.c[0].r;
			dg1 = t.c[1].g - t.c[0].g;
			db1 = t.c[1].b - t.c[0].b;
		}

		// (C.y - A.y > 0)
		if (t.vt[2].y - t.vt[0].y > 0) 
		{
			dx2=(t.vt[2].x - t.vt[0].x)/(t.vt[2].y - t.vt[0].y);
			dr2=(t.c[2].r - t.c[0].r)/(t.vt[2].y - t.vt[0].y);
			dg2=(t.c[2].g - t.c[0].g)/(t.vt[2].y - t.vt[0].y);
			db2=(t.c[2].b - t.c[0].b)/(t.vt[2].y - t.vt[0].y);
		} else 
			dx2=dr2=dg2=db2=0;

		// (C.y - B.y > 0)
		if (t.vt[2].y - t.vt[1].y > 0) 
		{
			dx3=(t.vt[2].x - t.vt[1].x)/(t.vt[2].y - t.vt[1].y);
			dr3=(t.c[2].r - t.c[1].r)/(t.vt[2].y - t.vt[1].y);
			dg3=(t.c[2].g - t.c[1].g)/(t.vt[2].y - t.vt[1].y);
			db3=(t.c[2].b - t.c[1].b)/(t.vt[2].y - t.vt[1].y);
		} else 
			dx3=dr3=dg3=db3=0;

		Vector3 S = t.vt[0];
		Vector3 E = t.vt[0];

		Color cs = t.c[0];
		Color ce = t.c[0];

		float dr, dg, db;
		Vector3 p;
		Color c;

		if(dx1 > dx2) 
		{
			for(; S.y <= t.vt[1].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;
				
				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx2; cs.r += dr2; cs.g += dg2; cs.b += db2;
				E.x += dx1; ce.r += dr1; ce.g += dg1; ce.b += db1;
			}

			E = t.vt[1]; 
			ce = t.c[1];

			for(; S.y <= t.vt[2].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx2; cs.r += dr2; cs.g += dg2; cs.b += db2;
				E.x += dx3; ce.r += dr3; ce.g += dg3; ce.b += db3;
			}

		}else
		{
			for(; S.y <= t.vt[1].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx1; cs.r += dr1; cs.g += dg1; cs.b += db1;
				E.x += dx2; ce.r += dr2; ce.g += dg2; ce.b += db2;
			}

			S = t.vt[1]; 
			cs = t.c[1];

			for(; S.y <= t.vt[2].y; S.y++, E.y++)
			{
				if(E.x - S.x > 0)
				{
					dr = (ce.r - cs.r)/(E.x - S.x);
					dg = (ce.g - cs.g)/(E.x - S.x);
					db = (ce.b - cs.b)/(E.x - S.x);
				}else
					dr=dg=db=0;

				p = S;
				c = cs;

				for(;p.x < E.x; p.x++)
				{
					putPixel(p, c);
					c.r += dr; c.g += dg; c.b += db;
				}
				S.x += dx3; cs.r += dr3; cs.g += dg3; cs.b += db3;
				E.x += dx2; ce.r += dr2; ce.g += dg2; ce.b += db2;
			}
		}
	}

//DRAW WIREFRAME
	if(state.style == Wireframe)
	{
		drawLine( t.c[0].r, t.c[0].g, t.c[0].b, t.vt[0].x, t.vt[0].y, t.c[1].r, t.c[1].g, t.c[1].b, t.vt[1].x, t.vt[1].y );
		drawLine( t.c[1].r, t.c[1].g, t.c[1].b, t.vt[1].x, t.vt[1].y, t.c[2].r, t.c[2].g, t.c[2].b, t.vt[2].x, t.vt[2].y );
		drawLine( t.c[0].r, t.c[0].g, t.c[0].b, t.vt[0].x, t.vt[0].y, t.c[2].r, t.c[2].g, t.c[2].b, t.vt[2].x, t.vt[2].y );		
	}

}

void Pipeline::bezierCurve(Vector3 *points, int order, int LOD)
{
	//DRAW BEZIER CURVE

	for(int i=0;i!=LOD;++i) 
	{

		// use the parametric time value 0 to 1
		float t = (float)i/(LOD);

		// pre-calculate 1.0f-t
		float it = 1.0f-t;
	
		// calculate blending functions
		float b0 = t*t*t;
		float b1 = 3*t*t*it;
		float b2 = 3*t*it*it;
		float b3 =  it*it*it;

		//Third order
		float x = b0*points[0].x +
				  b1*points[1].x +
				  b2*points[2].x +
				  b3*points[3].x ;

		float y = b0*points[0].y +
				  b1*points[1].y +
				  b2*points[2].y +
				  b3*points[3].y ;

		float z = b0*points[0].z +
				  b1*points[1].z +
				  b2*points[2].z +
				  b3*points[3].z ;
		
		// specify the point
		putPixel( Vector3(x,y,z), Color(1,0,1) );
	}

}