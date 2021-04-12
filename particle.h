#pragma once
#include "model.h"
class Particle
{
private:
	double mass, rad; //Mass and radius

public:

	Particle();

	vec3 pos; //Position radius-vector
	vec3 vel; //Velosity vector
	vec3 ax; //Axeleration vector
	double energy;

	double density;
	double pressure;

	double h = 3; //Smooth radius
	double gamma = 1;
	double A = 1; //Const dep on env

	double get_mass();
	double get_rad();

	void set_prop(double, double);

	void set_pos(double, double, double);
	void set_pos(vec3);

	void set_vel(double, double, double);
	void set_vel(vec3);

	void set_ax(double, double, double);
	void set_ax(vec3);

	void refresh();

};