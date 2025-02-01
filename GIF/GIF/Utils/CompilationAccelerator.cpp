#include "stdafx.h"

#include "CompilationAccelerator.h"



void CompilationAccelerator::insert_non_intersecting_curves(Arrangement_2& arr, const std::list<ESegment_2>& segmentsList)
{
	CGAL::insert_non_intersecting_curves(arr, segmentsList.begin(), segmentsList.end());
}


double CompilationAccelerator::to_double(const ARRNumberType& num)
{
	return CGAL::to_double(num);
}


bool CompilationAccelerator::assign(Arrangement_2::Vertex_const_handle& v, CGAL::Object& obj)
{
	return CGAL::assign(v, obj);
}


bool CompilationAccelerator::assign(Arrangement_2::Vertex_handle& v, CGAL::Object& obj)
{
	return CGAL::assign(v, obj);
}


bool CompilationAccelerator::assign(Arrangement_2::Halfedge_const_handle& h, CGAL::Object& obj)
{
	return CGAL::assign(h, obj);
}


bool CompilationAccelerator::assign(Arrangement_2::Halfedge_handle& h, CGAL::Object& obj)
{
	return CGAL::assign(h, obj);
}


bool CompilationAccelerator::assign(Arrangement_2::Face_const_handle& f, CGAL::Object& obj)
{
	return CGAL::assign(f, obj);
}


bool CompilationAccelerator::assign(Arrangement_2::Face_handle& f, CGAL::Object& obj)
{
	return CGAL::assign(f, obj);
}

bool CompilationAccelerator::isSimplePolygon(const std::vector<Complex>& polygon, bool bCheckOrientation)
{
	EPolygon_2 exactPolygon;
	int numVertices = polygon.size();

	for (int i = 0; i < numVertices; i++)
	{
		exactPolygon.push_back(EPoint_2(polygon[i].real(), polygon[i].imag()));
	}
	return isSimplePolygon(exactPolygon, bCheckOrientation);
}

bool CompilationAccelerator::isSimplePolygonInexact(const std::vector<Complex>& polygon, bool bCheckOrientation)
{
	typedef CGAL::Simple_cartesian<double> Kernel;
	CGAL::Polygon_2<Kernel> poly;
	int numVertices = polygon.size();

	for (int i = 0; i < numVertices; i++)
	{
		poly.push_back(Point_2(polygon[i].real(), polygon[i].imag()));
	}
	return isSimplePolygon(poly, bCheckOrientation);
}

//this function is a wrapper to triangulatePolygon that also uses CGAL to check if the polygon is simple before applying the "triangle" code for triangulating the polygon
bool CompilationAccelerator::triangulatePolygonSafe(const std::vector<std::complex<double> >& polygonPoints, std::vector<std::complex<double> >& meshVertices,
	std::vector<unsigned int>& triangleIndices, std::vector<std::pair<int, std::complex<double> > >& boundaryVertices,
	double maxTriangleArea, bool subsampleBoundaryEdges)
{
	bool isSimple = isSimplePolygon(polygonPoints);
	if (isSimple)
	{
		return triangulatePolygon(polygonPoints, meshVertices, triangleIndices, boundaryVertices, maxTriangleArea, subsampleBoundaryEdges);
	}
	else
	{
		return false;
	}
}


void CompilationAccelerator::zone(Arrangement_2& arr, const ESegment_2& seg, std::list<CGAL::Object>& intersectionList, Landmarks_pl& pointLocatorStructure)
{
	CGAL::zone(arr, seg, std::back_inserter(intersectionList), pointLocatorStructure);
}


void CompilationAccelerator::compute_intersection_points(const std::vector<ESegment_2>& segments, std::vector<EPoint_2>& intersectionPoints, bool report_endpoints)
{
	CGAL::compute_intersection_points(segments.begin(), segments.end(), back_inserter(intersectionPoints), report_endpoints);
}


/*

void CompilationAccelerator::overlay(const Arrangement_2& arr_in1, const Arrangement_2& arr_in2, Arrangement_2& arr_result, Overlay_traits& overlay_traits)
{
	CGAL::overlay(arr_in1, arr_in2, arr_result, overlay_traits);
}





*/
