#include "model.h"

Neighbour::Neighbour(size_t size_)
{
	this->resize(size_);
}

void Neighbour::init(std::vector<Particle>& particle)
{

	for (size_t i = 0; i < particle.size(); ++i)
	{
		for (size_t j = i+1; j < particle.size(); ++j)
		{

			//if (abs(particle[i].pos - particle[j].pos) < 2*particle[0].h) 
			//faster variant
			if ((particle[i].pos - particle[j].pos).abssq() < 4 * particle[i].h * particle[i].h)
			{
				web[i].push_back(j);
				web[j].push_back(i);
			}

		}
	}

}


std::vector<size_t> Neighbour::operator()(size_t n)
{
	return web[n];
}

void Neighbour::resize(size_t size_)
{
	size = size_;
	web.resize(size);
}

void Neighbour::fitsize(int x, int y)
{
	if (abs(x) > size)
		this->resize((size_t)abs(x) * 2 + 4);
	if (abs(y) > size)
		this->resize((size_t)abs(y) * 2 + 4);
}