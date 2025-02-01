#include "Utils\Utilities.h"

#include <ctime>
#include <limits>
#include <assert.h>

void swap(double& a, double& b)
{
	double temp = b;
	b = a;
	a = temp;
}

void swap(void** a, void** b)
{
	void* temp = *b;
	*b = *a;
	*a = temp;
}

int getRandomNumber()
{
	static bool firstTime = true;
	if(firstTime)
	{
		srand((unsigned)time(NULL));
		firstTime = false;
	}
	int randNumber = rand();
	return randNumber;
}


double getRandomNumber(double low, double high)
{
	static bool firstTime = true;
#define NORMALIZATION_RAND_FACTOR (3.0518509475997192297128208258309e-5) //(1.0 / RAND_MAX)
	if(firstTime)
	{
		srand((unsigned)time(NULL));
		firstTime = false;
	}
	double randNumber = (high - low)*(NORMALIZATION_RAND_FACTOR*rand()) + low;
	return randNumber;
}

double computePolygonArea(const std::vector<std::complex<double> >& polygon)
{
	int numVertices = polygon.size();
	double doubleArea = 0.0;

	std::complex<double> center(0.0, 0.0);
	for(int i = 0; i < numVertices; i++)
	{
		center += polygon[i];
	}
	for(int i = 0; i < numVertices; i++)
	{
		std::complex<double> v1 = polygon[i] - center;
		std::complex<double> v2 = polygon[i < numVertices - 1 ? i+1 : 0] - center;

		doubleArea += (v1.real()*v2.imag() - v1.imag()*v2.real());
	}
	return 0.5*doubleArea;
}