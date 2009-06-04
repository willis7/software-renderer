#pragma once
#include <math.h>

class Vector3
{
public:

	float x,y,z;

//Constructors

	//Default constructors leaves vector in
	//an intermediate state
	Vector3(){}

	//Copy constructor
	Vector3(const Vector3 &a) : x(a.x), y(a.y), z(a.z){}

	// Construct given three values
	Vector3(float nx, float ny, float nz) : x(nx), y(ny), z(nz) {}

//Standard object maintainance

	//Assignment
	Vector3 &operator = (const Vector3 &a)
	{
		x = a.x; y = a.y; z = a.z;
		return *this;
	}

	// Check for equality
	bool operator ==(const Vector3 &a) const {
		return x==a.x && y==a.y && z==a.z;
	}

	bool operator !=(const Vector3 &a) const 
	{
		return x!=a.x || y!=a.y || z!=a.z;
	}

// Vector operations

	// Vector operations
	// Set the vector to zero
	void zero() { x = y = z = 0.0f; }

	void set(float nx, float ny, float nz){x = nx; y = ny; z = nz;}

	// Unary minus returns the negative of the vector
	Vector3 operator -() const { return Vector3(-x,-y,-z); }

	// Binary + and – add and subtract vectors
	Vector3 operator +(const Vector3 &a) const {
		return Vector3(x + a.x, y + a.y, z + a.z);
	}

	Vector3 operator -(const Vector3 &a) const {
		return Vector3(x - a.x, y - a.y, z - a.z);
	}

	// Multiplication and division by scalar
	Vector3 operator *(float a) const {
		return Vector3(x*a, y*a, z*a);
	}

	Vector3 operator /(float a) const {
		float oneOverA = 1.0f / a;
		return Vector3(x*oneOverA, y*oneOverA, z*oneOverA);
	}

	Vector3 &operator +=(const Vector3 &a) {
		x += a.x; y += a.y; z += a.z;
		return *this;
	}

	Vector3 &operator -=(const Vector3 &a) {
		x -= a.x; y -= a.y; z -= a.z;
		return *this;
	}

	Vector3 &operator*=(float a) {
		x *= a; y *= a; z *= a;
		return *this;
	}

	Vector3 &operator /=(float a) {
		float oneOverA = 1.0f / a;
		x *= oneOverA; y *= oneOverA; z *= oneOverA;
		return *this;
	}

	// Normalize the vector
	void normalize() {
		float magSq = x*x + y*y + z*z;
		if (magSq > 0.0f) { // check for divide-by-zero
		float oneOverMag = 1.0f / sqrt(magSq);
		x *= oneOverMag;
		y *= oneOverMag;
		z *= oneOverMag;
		}
	}

	// Vector dot product. Overloaded 
	// multiplication symbol
	float operator *(const Vector3 &a) const {
		return x*a.x + y*a.y + z*a.z;
	}

};

// Compute the length of a vector
inline float vecLength(const Vector3 &a) {
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

// Compute the cross product of two vectors
inline Vector3 crossProd(const Vector3 &a, const Vector3 &b) {
	return Vector3(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	);
}

// Scalar on the left multiplication, for symmetry
inline Vector3 operator *(float k, const Vector3 &v) {
	return Vector3(k*v.x, k*v.y, k*v.z);
}

// Compute the distance between two points
inline float distance(const Vector3 &a, const Vector3 &b) {
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float dz = a.z - b.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

// Global zero vector constant
extern const Vector3 kZeroVector;


// Represents a vertex in homogenous coordinates
class Vector4
{
public:

	float x,y,z,w;

//Constructors

	//Default constructors leaves vector in
	//an intermediate state
	Vector4(){}

	//Copy constructor
	Vector4(const Vector4 &a) : x(a.x), y(a.y), z(a.z){}

	// Construct given three values
	Vector4(float nx, float ny, float nz) : x(nx), y(ny), z(nz) { w = 1;}

	Vector4(float nx, float ny, float nz, float nw) : x(nx), y(ny), z(nz), w(nw){}

//Standard object maintainance

	//Assignment
	Vector4 &operator = (const Vector4 &a)
	{
		x = a.x; y = a.y; z = a.z;
		return *this;
	}

	// Check for equality
	bool operator ==(const Vector4 &a) const {
		return x==a.x && y==a.y && z==a.z;
	}

	bool operator !=(const Vector4 &a) const 
	{
		return x!=a.x || y!=a.y || z!=a.z;
	}

// Vector operations

	// Vector operations
	// Set the vector to zero
	void zero() { x = y = z = w = 0.0f; }

	// Unary minus returns the negative of the vector
	Vector4 operator -() const { return Vector4(-x,-y,-z); }

	// Binary + and – add and subtract vectors
	Vector4 operator +(const Vector4 &a) const {
		return Vector4(x + a.x, y + a.y, z + a.z);
	}

	Vector4 operator -(const Vector4 &a) const {
		return Vector4(x - a.x, y - a.y, z - a.z);
	}

	// Multiplication and division by scalar
	Vector4 operator *(float a) const {
		return Vector4(x*a, y*a, z*a);
	}

	Vector4 operator /(float a) const {
		float oneOverA = 1.0f / a;
		return Vector4(x*oneOverA, y*oneOverA, z*oneOverA);
	}

	Vector4 &operator +=(const Vector4 &a) {
		x += a.x; y += a.y; z += a.z;
		return *this;
	}

	Vector4 &operator -=(const Vector4 &a) {
		x -= a.x; y -= a.y; z -= a.z;
		return *this;
	}

	Vector4 &operator*=(float a) {
		x *= a; y *= a; z *= a;
		return *this;
	}

	Vector4 &operator /=(float a) {
		float oneOverA = 1.0f / a;
		x *= oneOverA; y *= oneOverA; z *= oneOverA;
		return *this;
	}

	// Normalize the vector
	void normalize() {
		float magSq = x*x + y*y + z*z;
		if (magSq > 0.0f) { // check for divide-by-zero
		float oneOverMag = 1.0f / sqrt(magSq);
		x *= oneOverMag;
		y *= oneOverMag;
		z *= oneOverMag;
		}
	}

	// Vector dot product. Overloaded 
	// multiplication symbol
	float operator *(const Vector4 &a) const {
		return x*a.x + y*a.y + z*a.z;
	}

};

// Compute the length of a vector
inline float vecLength(const Vector4 &a) {
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

// Compute the cross product of two vectors
inline Vector4 crossProd(const Vector4 &a, const Vector4 &b) {
	return Vector4(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	);
}

// Scalar on the left multiplication, for symmetry
inline Vector4 operator *(float k, const Vector4 &v) {
	return Vector4(k*v.x, k*v.y, k*v.z);
}

// Compute the distance between two points
inline float distance(const Vector4 &a, const Vector4 &b) {
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float dz = a.z - b.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

