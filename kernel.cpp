//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#include "model.h"



int main()
{

	const int N = 6 * 18 + 4 * 6;
	double dt = 0.0025/4;
	
	double t = 6;

	Kernel kernel(2);
	std::vector<Particle> particle(N);


	int c = 0;
	for(double i = 0; i<6; ++i)
		for (double j = 0; j < 18; ++j)
		{
			if((int)i % 2 == 0)
				particle[c].set_pos(i*0.86, j + 0.5, 0);
			else
				particle[c].set_pos(i * 0.86, j, 0);
			//particle[c].set_vel((rand() % 2 - 1) * 1, (rand() % 2 - 1) * 1, 0);
			particle[c].set_prop(100, 1);
			++c;
			system("cls");
			std::cout << "Getenrating particles: " << c << "/" << N << std::endl;
		}

	
	for (double i = 0; i < 6; ++i)
		for (double j = 0; j < 4; ++j)
		{
			particle[c].set_pos(i + 10, j + 6, 0);
			particle[c].set_vel( -10, 0, 0);
			particle[c].set_prop(25, 1);
			++c;
			system("cls");
			std::cout << "Getenrating particles: " << c << "/" << N << std::endl;
		}
	
	
	std::ofstream file("test.txt");
	//std::ofstream energy("data_6_en");

	size_t i = 0;
	//double energy = 0;
	unsigned int start_time = clock();
	while (i < t / dt)
	{

		//advanced_scheme(particle, kernel, dt, i);
		eiler_scheme(particle, kernel, dt, i);

		if (i % (int)(t * 0.1 / dt) == 0)
		{
			system("cls");
			std::cout << "Progress: " << i * dt * 100 / t << " % "<<std::endl
					  << "interior time = " << i*dt << " s";
		}

		for (auto& p : particle)
			file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " ";
		file << std::endl;

		++i;
	}
	unsigned int end_time = clock();
	std::cout << std::endl << (end_time - start_time) / i;

	file.close();

}
