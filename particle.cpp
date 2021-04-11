#include "model.h"

Particle::Particle()
{
	set_prop(1, 1);
	pos.set(0, 0, 0);
	vel.set(0, 0, 0);
	ax.set(0, 0, 0);
}

double Particle::get_mass()
{
	return mass;
}

double Particle::get_rad()
{
	return rad;
}

void Particle::set_prop(double m, double R)
{
	mass = m;
	rad = R;
}

void Particle::set_pos(double X, double Y, double Z)
{
	pos.set(X, Y, Z);
}
void Particle::set_pos(vec3 other)
{
	pos.set(other);
}

void Particle::set_vel(double X, double Y, double Z)
{
	vel.set(X, Y, Z);
}
void Particle::set_vel(vec3 other)
{
	vel.set(other);
}

void Particle::set_ax(double X, double Y, double Z)
{
	ax.set(X, Y, Z);
}
void Particle::set_ax(vec3 other)
{
	ax.set(other);
}


