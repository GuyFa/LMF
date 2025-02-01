
#pragma once

namespace CGUT
{
	double length(const Vector_3& v);
	void normalize(Vector_3& v);
	double exactFaceOrientation(Facet_const_handle& face);
	bool isPointInsideFaceUV(Facet_const_handle& face, const std::complex<double>& p);
	//Point_3 mapPointInsideTriangle(Facet_const_handle& face, const Complex& uv);
	int countFoldovers(const Mesh& mesh, bool printReport = true);
	void point3ToComplex(const Point_3& p, std::complex<double>& c);
	const std::complex<double> point3ToComplex(const Point_3& p);
	void complexToPoint3(const std::complex<double>& c, Point_3& p);
	const Point_3 complexToPoint3(const std::complex<double>& c);
}


