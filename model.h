#pragma once
#pragma once
#include "Vmath.h" 
#include "SPH.h"
#include "particle.h"

#include <omp.h>
#include <functional>

//Functor to find neighbours for the particle
class Neighbour
{
public:
	std::vector<std::vector<size_t>> web;
	size_t size;

	Neighbour(size_t);

	void init(std::vector<Particle>&);
	void resize(size_t);
	void fitsize(int, int);

	//Returns vector with neighbours' numbers
	std::vector<size_t> operator()(size_t n);

};


//Computing schemes
void eiler_scheme(std::vector<Particle>&, Kernel&, double, size_t, std::function<bool(vec3)>);
void advanced_scheme(std::vector<Particle>&, Kernel&, double, size_t);

bool nei(Particle&, Particle&);


// particle's axeleration, dv/dt 
vec3 ax(size_t, std::vector<Particle>&, Kernel&, Neighbour&);
vec3 ax_inv_Eu(size_t, std::vector<Particle>&, Kernel&, Neighbour&);
double energy(std::vector<Particle>&); 


//density of smooth particle (particle number, vector with particles, kernel)
double dens(size_t, std::vector<Particle>&, Kernel&);
double dens(size_t n, std::vector<Particle>&, Kernel&, Neighbour&);

//pressure (-||-)
double press(size_t, double, std::vector<Particle>&, Kernel&, Neighbour&);
//adiabatic exponent (-||-)
double gamma(size_t, double, std::vector<Particle>&, Kernel&); 
//adaptive smooth radius  (-||-)
void adapt_h(size_t, double, std::vector<Particle>&, Kernel&); 

void refresh(std::vector<Particle>&, Kernel&, Neighbour&);

//rotor and divergence of velocity's field (-||-)
double  divVel(size_t, std::vector<Particle>&, Kernel&);
vec3	rotVel(size_t, std::vector<Particle>&, Kernel&);

//Lennard_Johns potencial field
double U_lj(vec3, vec3);




