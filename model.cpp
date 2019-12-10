#include "Vmath.h"
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




double dens(size_t n, std::vector<Particle>& particle, Kernel& kernel)
{
	double dens = 0;
	for (Particle& p : particle)
		if (nei(particle[n], p))
			dens += p.get_mass()* kernel.W(particle[n].pos, p.pos, particle[n].h);

	return dens;
}
double dens(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	double dens = 0;
	for (size_t i : neighbour(particle[n]))
		if (nei(particle[n], particle[i]))
			dens += particle[i].get_mass() * kernel.W(particle[n].pos, particle[i].pos, particle[n].h);

	return dens;
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
	for (size_t i : nei(particle[n]))
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
	for (size_t i : nei(particle[n]))
		rotV += (particle[i].vel - particle[n].vel) /
		kernel.gradW(particle[n].pos, particle[i].pos, particle[n].h) *
		particle[i].get_mass();

	return rotV / dens(n, particle, kernel);
}



void eiler_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration)
{
	Neighbour neighbour(10);
	neighbour.init(particle);
	size_t n = 0;
	for (Particle& p : particle)
	{	
		p.ax = ax(n, particle, kernel, neighbour);
		p.vel += p.ax * dt;
		p.pos += p.vel * dt;

		++n;
	}
}

void advanced_scheme(std::vector<Particle>& particle, Kernel& kernel, double dt, size_t iteration)
{	
	size_t n = 0;
	std::vector<vec3> v_0;
	std::vector<vec3> pos_0;
	std::vector<vec3> pos_part;

	Neighbour neighbour(10);
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

	neighbour.init(particle);
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

vec3 ax(size_t n, std::vector<Particle>& particle, Kernel& kernel, Neighbour& neighbour)
{
	vec3 ax (0, 0, 0);
	
	for (size_t i : neighbour(particle[n]))
	{
		if (i != n && nei(particle[n], particle[i]))
			ax += kernel.gradW(particle[n].pos, particle[i].pos, 1) *
			particle[i].get_mass() * U_lj(particle[n].pos, particle[i].pos) / dens(n, particle, kernel, neighbour);
		else continue;
	}
	return ax * (-1) / particle[n].get_mass();
}

double U_lj(vec3 r1, vec3 r2)
{	
	double s = 0.891;
	double e = 1;
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


Neighbour::Neighbour(size_t size_)
{	
	this->resize(size_);
}

void Neighbour::init(std::vector<Particle>& particle)
{
	for (auto& w1 : web)
		for (auto& w2 : w1)
			w2.resize(0);
	
	int x;
	int y;
	for (size_t i = 0; i < particle.size(); ++i)
	{
		x = (int)(particle[i].pos.x / particle[i].h);
		y = (int)(particle[i].pos.y / particle[i].h); 
		fitsize(x, y);

		web[x + size][y + size].push_back(i);
	}

}

std::vector<size_t> Neighbour::operator()(Particle& p)
{
	std::vector<size_t> res(0);
	int x = (int)(p.pos.x / p.h);
	int y = (int)(p.pos.y / p.h);
	fitsize(x, y);

	x += size;
	y += size;
	for (int i = -2; i < 3; ++i)
		for (int j = -2; j < 3; ++j)
			res.insert(res.end(), web[(size_t)(x + i)][y + j].begin(), web[x + i][y + j].end());
	
	return res;
}

void Neighbour::resize(size_t size_)
{
	size = size_;
	web.resize(size * 2);
	for (auto& num : web)
		num.resize(size * 2);
}

void Neighbour::fitsize(int x, int y)
{
	if (abs(x) > size)
		this->resize((size_t)abs(x) * 2 + 4);
	if (abs(y) > size)
		this->resize((size_t)abs(y) * 2 + 4);
}