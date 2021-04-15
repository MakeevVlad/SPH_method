#include "model.h"



double dens(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double dens = 0;
	for (Particle& p : particle)
		dens += p.get_mass()* kernel.W(particle[n].pos, p.pos, particle[n].h);

	return dens;
}
double dens(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	double dens = 0;
#pragma omp parallel for reduction(+:dens)
	for (size_t i : neighbour(n))
		dens += particle[i].get_mass() * kernel.W(particle[n].pos, particle[i].pos, particle[n].h);

	return dens;
}
double press(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	return particle[n].A * pow(particle[n].density, particle[n].gamma);
}

void adapt_h(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel)
{
	particle[n].h = particle[n].h * (1 + divVel(n, particle, kernel)/kernel.D);
}

double divVel(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double divV = 0;
	for (Particle& p : particle)
		if (nei(particle[n], p))
			divV += (p.vel - particle[n].vel) * 
					kernel.gradW(particle[n].pos, p.pos, particle[n].h) * 
					p.get_mass();

	return divV / dens(n, particle, kernel);

}
double divVel(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& nei)
{
	double divV = 0;
	for (size_t i : nei(n))
		divV += (particle[i].vel - particle[n].vel) *
		kernel.gradW(particle[n].pos, particle[i].pos, particle[n].h) *
		particle[i].get_mass();

	return divV / dens(n, particle, kernel);

}

vec3 rotVel(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel)
{
	vec3 rotV = 0;
	for (Particle& p : particle)
		if (nei(particle[n], p))
			rotV += (p.vel - particle[n].vel) /
			kernel.gradW(particle[n].pos, p.pos, particle[n].h) *
			p.get_mass();

	return rotV / dens(n, particle, kernel);
}vec3 rotVel(size_t n, double dt, std::vector<Particle>& particle, Kernel& kernel, Neighbour& nei)
{
	vec3 rotV = 0;
	for (size_t i : nei(n))
		rotV += (particle[i].vel - particle[n].vel) /
		kernel.gradW(particle[n].pos, particle[i].pos, particle[n].h) *
		particle[i].get_mass();

	return rotV / dens(n, particle, kernel);
}

void refresh(std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	for (size_t i = 0; i < particle.size(); ++i)
	{
		particle[i].density = dens(i, particle, kernel, neighbour);
		particle[i].pressure = press(i, particle, kernel, neighbour);
	}
}



vec3 ax(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	

	vec3 ax (0, 0, 0);
	#pragma omp parallel for reduction(+:ax)
	for (size_t i : neighbour(n))
	{
		if (i != n )
			ax += kernel.gradW(particle[n].pos, particle[i].pos, particle[n].h) *
			particle[i].get_mass() * U_lj(particle[n].pos, particle[i].pos) / particle[i].density;
	}
	return ax * (-1) / particle[n].get_mass();
}

vec3 ax_inv_Eu(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	vec3 ax(0, 0, 0);
	vec3 g(0, -5, 0);

#pragma omp parallel for reduction(+:ax)
	for (size_t i : neighbour(n))
	{
		if (i != n)
		{
			ax += kernel.gradW(particle[n].pos, particle[i].pos, particle[n].h) *(particle[n].pressure / particle[n].density / particle[n].density +
				 particle[i].pressure / particle[i].density / particle[i].density) * particle[i].get_mass();
		}
	}
	return ax * (-1) / particle[n].get_mass() + g;
}

double U_lj(vec3 r1, vec3 r2)
{	
	double s = 0.89;
	double e = 4;
	double r = abs(r1 - r2);
	double s_r = pow(s / r, 6);

	return 4 * e * (s_r * s_r - s_r);
}

bool nei(Particle& p0, Particle& p1)
{
	if (abs(p1.pos.x - p0.pos.x) < 2 * p0.h
		&& abs(p1.pos.y - p0.pos.y) < 2 * p0.h
		&& abs(p1.pos.z - p0.pos.z) < 2 * p0.h)
		return 1;
	else
		return 0;
}
