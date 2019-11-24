#include "SPH.h"

Kernel::Kernel(size_t D_)
{
	D = D_;
	if (D == 1) C = 0.6667; // = 2/3
	else if (D == 2) C = 0.4547; // = 10/(7pi)
	else if (D == 3) C = 0.3183; // = 1/pi

}

double h_value()
{
	return 1;
}

double Kernel::W(vec3 r1, vec3 r2, double h)
{
	double k = abs(r1 - r2) / h ;
	if (k <= 1)
	{
		return C * (1 - 1.5 * pow(k, 2) + 0.75 * pow(k, 3)) / pow(h, D);
	}
	else if (k > 1 && k <= 2)
	{
		return C * 0.25 * pow(2 - k, 3) / pow(h, D);
	}
	else return 0;

}

vec3 Kernel::gradW(vec3 r1, vec3 r2, double h)
{
	double k = abs(r1 - r2) / h;
	if (k <= 1)
	{
		return (r1 - r2) * C * (-3*k + 2.25*pow(k, 2))/ (pow(h, D) * abs(r1 - r2));
	}
	else if (k > 1 && k <= 2)
	{
		return (r1 - r2) * C * 0.75 * pow(2 - k, 2) / (pow(h, D) * abs(r1 - r2));
	}
	else return vec3(0, 0, 0);
}