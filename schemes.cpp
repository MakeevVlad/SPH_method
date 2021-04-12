#include "model.h"

void eiler_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration, std::function<bool(vec3)> bounds)
{
	Neighbour neighbour(particle.size());
	neighbour.init(particle);

	refresh(particle, kernel, neighbour);

	size_t n = 0;
	#pragma omp parallel for
	for (Particle& p : particle)
	{
		p.ax = ax_inv_Eu(n, particle, kernel, neighbour);
		++n;
	}
	for (Particle& p : particle)
	{
		p.vel += p.ax * dt;
		p.pos += p.vel * dt;
		if (bounds(p.pos))
		{
			p.vel = p.vel * (-1);
			p.pos = p.vel * dt * (-1);
		}
	}
}

void advanced_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration)
{
	size_t n = 0;
	std::vector<vec3> v_0;
	std::vector<vec3> pos_0;
	std::vector<vec3> pos_part;

	Neighbour neighbour(particle.size());
	neighbour.init(particle);

	std::vector<vec3> axeleration(0);
	for (Particle& p : particle)
	{
		p.ax = ax(n, particle, kernel, neighbour);
		p.vel = p.vel + p.ax * dt * 0.5;

		pos_0.push_back(p.pos);
		v_0.push_back(p.vel);

		pos_part.push_back(p.pos + (p.vel + v_0[n]) * dt * 0.25);

		++n;
	}
	for (n = 0; n < particle.size(); ++n)
	{
		particle[n].pos = pos_part[n];
	}

	for (n = 0; n < particle.size(); ++n)
	{
		particle[n].ax = pos_part[n];
		particle[n].vel = v_0[n] + particle[n].ax * dt;
	}
	for (n = 0; n < particle.size(); ++n)
	{
		particle[n].pos = pos_0[n] + (v_0[n] + particle[n].vel) * 0.5 * dt;
	}

}