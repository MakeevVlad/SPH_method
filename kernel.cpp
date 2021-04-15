//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include "model.h"

bool bounds(vec3 pos)
{
	return ((pos.x > 11) || (pos.x < -1) || (pos.y > 11) || (pos.y < -1));
}


int main()
{

	const int N = 400;
	double k = 10 / sqrt(N);
	double dt = 0.006;
	
	double t = 20;

	Kernel kernel(2);
	std::vector<Particle> particle(N);


	int c = 0;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-5, 5);
	for (double i = 0; i < sqrt(N); ++i)
		for (double j = 0; j < sqrt(N); ++j)
		{
			particle[c].set_pos(i*k, j*k, 0);
			particle[c].set_vel( dis(gen), dis(gen), 0);
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
		eiler_scheme(particle, kernel, dt, i, bounds);

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
	std::cout << std::endl << "Full time: " << (end_time - start_time);
	std::cout << std::endl << "Mean iteration time: " << (end_time - start_time) / i;

	file.close();

}
