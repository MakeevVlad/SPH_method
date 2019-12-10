//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Vmath.h"
#include "Model.h"



int main()
{

	const int N = 400;
	double dt = 0.0025;
	
	double t = 1;

	Kernel kernel(2);
	std::vector<Particle> particle(N);

	double x = -10;
	double y = -10;

	size_t c = 0;
	for (size_t i = 0; i < particle.size(); ++i)
	{

		if (i % 20 == 0 && i != 0)
		{
			x = -10; 
			y += 1;
			++c;
		}
		particle[i].set_pos(x, y, 0);
		//particle[i].set_vel((rand() % 10 - 5) * 0.01, (rand() % 10 - 5) * 0.01, (rand() % 10 - 5) * 0.01);
		particle[i].set_prop(100, 1);
		x += 1;

		system("cls");
		std::cout << "Getenrating particles: " << i << "/" << N <<std::endl;
	}
	particle.resize(N + 4);
	particle[N].set_pos(-11, 0, 0);
	particle[N+1].set_pos(-12, 0, 0);
	particle[N+2].set_pos(-11, -1, 0);
	particle[N+3].set_pos(-12, -1, 0);
	particle[N].set_vel(11, 0, 0);
	particle[N + 1].set_vel(11, 0, 0);
	particle[N + 2].set_vel(11, 0, 0);
	particle[N + 3].set_vel(11, 0, 0);
	particle[N].set_prop(500, 1);
	particle[N + 1].set_prop(500, 1);
	particle[N + 2].set_prop(500, 1);
	particle[N + 3].set_prop(500, 1);
	
	std::ofstream file("data_19");
	//std::ofstream energy("data_6_en");

	size_t i = 0;
	//double energy = 0;
	while (i < t / dt)
	{

		advanced_scheme(particle, kernel, dt, i);


		system("cls");
		std::cout << "Progress: " << i * dt * 100 / t << " % "<<std::endl
				  << "interior time = " << i*dt << " s";
		

		for (auto& p : particle)
			file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " ";
		file << std::endl;

		++i;
	}

	file.close();

}
