#pragma once
#pragma once
#include "Vmath.h"
#include "SPH.h"




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

	
	double h = 2; //Smooth radius
	
	double get_mass();
	double get_rad();

	void set_prop(double, double);

	void set_pos(double, double, double);
	void set_pos(vec3);

	void set_vel(double, double, double);
	void set_vel(vec3);

	void set_ax(double, double, double);
	void set_ax(vec3);

};


//Functor to find neighbours for the particle
class Neighbour
{
public:
	std::vector<std::vector<std::vector<size_t>>> web;
	size_t size;

	Neighbour(size_t);

	void init(std::vector<Particle>&);
	void resize(size_t);
	void fitsize(int, int);

	//Returns vector with neighbours' numbers
	std::vector<size_t> operator()(Particle&);

};


//Computing schemes
void eiler_scheme(std::vector<Particle>&, Kernel&, double, size_t);
void advanced_scheme(std::vector<Particle>&, Kernel&, double, size_t);

bool nei(Particle&, Particle&);

// particle's axeleration, dv/dt 
vec3 ax(size_t, std::vector<Particle>&, Kernel&, Neighbour&);
double energy(std::vector<Particle>&); 


//density of smooth particle (particle number, vector with particles, kernel)
double dens(size_t, std::vector<Particle>&, Kernel&);
//pressure (-||-)
double pres(size_t, double, std::vector<Particle>&, Kernel&); 
//adiabatic exponent (-||-)
double gamma(size_t, double, std::vector<Particle>&, Kernel&); 
//adaptive smooth radius  (-||-)
void adapt_h(size_t, double, std::vector<Particle>&, Kernel&); 


//rotor and divergence of velocity's field (-||-)
double  divVel(size_t, std::vector<Particle>&, Kernel&);
vec3	rotVel(size_t, std::vector<Particle>&, Kernel&);

//Lennard_Johns potencial field
double U_lj(vec3, vec3);


