#pragma once
#include <iostream>
#include "Vector3.h"
using namespace std;

class Matrix4
{

public:

	float mat4[16];	

//Constructors

	//Default constructors sets to identity
	Matrix4(void);

	//Copy constructor
	Matrix4(const Matrix4 &matrix);

	//Destructor
	~Matrix4(void){}

//Matrix operations

	//Set to identity & check is identity
	void identity();
	bool isIdentity();

	// Assignment operator
	Matrix4& operator=(const Matrix4 &m);

	// Access index element
	float& operator[](const int index);

	//Matrix multiplication
	Matrix4 operator*(const Matrix4 &m)const;

	//Matrix transformation of 3D vector
	Vector3 operator*(const Vector3 &v) const;

	//Inverse matrix
	Matrix4& affineInverse();

//Object maintainance

	//Equality operators
	bool operator == (const Matrix4 &m) const;
	bool operator != (const Matrix4 &m) const;

//LookUp tables

	//flag to test for instantiation
	static bool tableCreated;

	//storage for look-up tables
	static float SIN_LUT[360];
	static float COS_LUT[360];

	//Fills the LUT arrays for later use
	void setupLUT();
	
	static float IDENTITY[16];

	
};
