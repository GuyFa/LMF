#pragma once

#include <complex>
#include <vector>


void swap(double& a, double& b);
void swap(void** a, void** b);
int getRandomNumber();
double getRandomNumber(double low, double high);
double computePolygonArea(const std::vector<std::complex<double> >& polygon);
//maps a scalar value to RGB color space based on Matlab's color jet. the mapping is linear.
//values equal or smaller than minValue are mapped to the blue color and values equal or larger than maxValue are mapped to the red color.
//values at the center of the range are mapped to the green color.
//see the jet colors here: http://www.mathworks.com/help/matlab/ref/colormap.html
void mapColor(double scalarValue, float& R, float& G, float& B, double minValue, double maxValue);



template<class T>
T P2(T x)
{
	return x*x;
}

template<class T>
T P3(T x)
{
	return x*x*x;
}