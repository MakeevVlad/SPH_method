#include "Vmath.h"

vec3::vec3(double* other) 
{ 
	set(other[0], other[1], other[2]); 
}

vec3::vec3(double X, double Y, double Z) 
{ 
	set(X, Y, Z); 
}

vec3::vec3() 
{ 
	set(0, 0, 0); 
}


void vec3::set(double X, double Y, double Z)
{
	x = X;
	y = Y;
	z = Z;
}
void vec3::set(vec3 other)
{
	set(other.x, other.y, other.z);
}

vec3 vec3::operator+ (vec3 other)
{
	vec3 res;
	res.x = x + other.x;
	res.y = y + other.y;
	res.z = z + other.z;
	return res;
}
vec3 vec3::operator+= (vec3& other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	return *this;
}
vec3 vec3::operator- (vec3 other)
{
	vec3 res;
	res.x = x - other.x;
	res.y = y - other.y;
	res.z = z - other.z;
	return res;
}
vec3 vec3::operator-= (vec3& other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	return *this;
}
vec3 vec3::operator*(double c)
{
	vec3 res(x * c, y * c, z * c);
	return res;
}
vec3 vec3::operator/(double c)
{
	vec3 res(x / c, y / c, z / c);
	return res;
}
double vec3::operator* (vec3 other)
{
	return x * other.x + y * other.y + z * other.z;
}
vec3 vec3::operator/(vec3 other)
{
	return vec3( y * other.z - z * other.y,
				 z * other.x - x * other.z,
				 x * other.y - y * other.x);
}

vec3& vec3::operator= (const vec3& other)
{
	x = other.x;
	y = other.y;
	z = other.z;
	return *this;
}
vec3 vec3::operator= (double* other)
{
	vec3 res(other);
	return res;
}

double vec3::abs()
{
	return sqrt(x * x + y * y + z * z);
}
double vec3::abssq()
{
	return x * x + y * y + z * z;
}
double vec3::proj(vec3 other)
{
	return this->abs() * cos(*this, other) ;
}



double abs(vec3 data)
{
	return sqrt(data.x * data.x + data.y * data.y + data.z * data.z);
}
double cos(vec3 v1, vec3 v2)
{
	return v1 * v2 / (v1.abs() * v2.abs());
}