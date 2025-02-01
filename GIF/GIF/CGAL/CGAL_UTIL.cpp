#include "Utils/stdafx.h"

#include "Utils/STL_Macros.h"
#include "CGAL_Mesh.h"
#include "CGAL_Macros.h"
#include "CGAL_UTIL.h"
#include "Utils/triangulation.h"

namespace CGUT
{

	typedef CGAL::Simple_cartesian<double>::Vector_3 Vector_3;
	typedef CGAL::Simple_cartesian<double>::Point_3  Point_3;
	typedef std::complex<double> Complex;

}


double CGUT::length(const Vector_3& v)
{
	return sqrt(v.squared_length());
}

void CGUT::normalize(Vector_3& v)
{
	double lensq = v.squared_length();

	if(lensq != 0.0)
	{
		v = v/sqrt(lensq);
	}
}


void CGUT::point3ToComplex(const Point_3& p, std::complex<double>& c)
{
	c.real(p.x());
	c.imag(p.y());
}


const std::complex<double> CGUT::point3ToComplex(const Point_3& p)
{
	std::complex<double> c(p.x(), p.y());
	return c;
}


void CGUT::complexToPoint3(const std::complex<double>& c, Point_3& p)
{
	p = Point_3(real(c), imag(c), 0.0);
}


const Point_3 CGUT::complexToPoint3(const std::complex<double>& c)
{
	Point_3 p(real(c), imag(c), 0.0);
	return p;
}

double CGUT::exactFaceOrientation(Facet_const_handle& face)
{
	Vertex_const_handle p[3];
	face->getVertices(p);

	Complex p1(p[0]->uv().x(), p[0]->uv().y());
	Complex p2(p[1]->uv().x(), p[1]->uv().y());
	Complex p3(p[2]->uv().x(), p[2]->uv().y());

	double cross = CounterClockWise(p1, p2, p3);

	return cross;
}

bool CGUT::isPointInsideFaceUV(Facet_const_handle& face, const std::complex<double>& uv)
{
	Vertex_const_handle v[3];
	face->getVertices(v);

	std::complex<double> uv1 = v[0]->uvC();
	std::complex<double> uv2 = v[1]->uvC();
	std::complex<double> uv3 = v[2]->uvC();

	bool isInside = CounterClockWise(uv1, uv2, uv) >= 0 && CounterClockWise(uv2, uv3, uv) >= 0 && CounterClockWise(uv3, uv1, uv) >= 0;

	return isInside;
}
/*

//assuming uv is inside the triangle in the uv domain. Find the corresponding point inside the triangle in 3D domain
Point_3 CGUT::mapPointInsideTriangle(Facet_const_handle& face, const Complex& uv)
{
	Vertex_const_handle v[3];
	face->getVertices(v);

	Complex uv1 = v[0]->uvC();
	Complex uv2 = v[1]->uvC();
	Complex uv3 = v[2]->uvC();



	double w1 = ((uv2.imag()-uv3.imag())*(uv.real()-uv3.real()) + (uv3.real()-uv2.real())*(uv.imag()-uv3.imag()))/((uv2.imag()-uv3.imag())*(uv1.real()-uv3.real()) + (uv3.real()-uv2.real())*(uv1.imag()-uv3.imag()));
	double w2 = ((uv3.imag()-uv1.imag())*(uv.real()-uv3.real()) + (uv1.real()-uv3.real())*(uv.imag()-uv3.imag()))/((uv3.imag()-uv1.imag())*(uv2.real()-uv3.real()) + (uv1.real()-uv3.real())*(uv2.imag()-uv3.imag()));
	double w3 = 1.0 - w1 - w2;

	Point_3 mappedPoint = CGAL::ORIGIN + (w1*(v[0]->point() - CGAL::ORIGIN) + w2*(v[1]->point() - CGAL::ORIGIN) + w3*(v[2]->point() - CGAL::ORIGIN));

	return mappedPoint;
}
*/

int CGUT::countFoldovers(const Mesh& mesh, bool printReport)
{
	assert(mesh.is_pure_triangle());
	int numFlipped = 0;
	int numGood = 0;

	for_each_const_facet(face, mesh)
	{
		double cross = exactFaceOrientation(face);

		if(cross <= 0.0)
		{
			numFlipped++;
		}
		else
		{
			numGood++;
		}
	}
	int numTotal = numFlipped + numGood;

	if(printReport)
	{
		cout << "There are: " << numFlipped << " flipped triangles out of " << numTotal << " in total" << endl;
	}

	return numFlipped;
}

