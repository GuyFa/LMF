#pragma once

#include "CGAL\CGAL_Mesh.h"
#include "CGAL\CGAL_Macros.h"
#include "Utils/STL_Macros.h"




class CompilationAccelerator
{
public:

	CompilationAccelerator() {};
	static void insert_non_intersecting_curves(Arrangement_2& arr, const std::list<ESegment_2>& segmentsList);
	// 	static void overlay(const Arrangement_2& arr_in1, const Arrangement_2& arr_in2, Arrangement_2& arr_result, Overlay_traits& overlay_traits);
	static void zone(Arrangement_2& arr, const ESegment_2& seg, std::list<CGAL::Object>& intersectionList, Landmarks_pl& pointLocatorStructure);
	static void compute_intersection_points(const std::vector<ESegment_2>& segments, std::vector<EPoint_2>& intersectionPoints, bool report_endpoints);
	static double to_double(const ARRNumberType& num);
	static bool assign(Arrangement_2::Vertex_const_handle& v, CGAL::Object& obj);
	static bool assign(Arrangement_2::Vertex_handle& v, CGAL::Object& obj);
	static bool assign(Arrangement_2::Halfedge_const_handle& h, CGAL::Object& obj);
	static bool assign(Arrangement_2::Halfedge_handle& h, CGAL::Object& obj);
	static bool assign(Arrangement_2::Face_const_handle& f, CGAL::Object& obj);
	static bool assign(Arrangement_2::Face_handle& f, CGAL::Object& obj);

	template<class K>
	static bool isSimplePolygon(const CGAL::Polygon_2<K>& polygon, bool bCheckOrientation = true)
	{
		if (bCheckOrientation)
			return (polygon.is_simple() && (polygon.orientation() == CGAL::COUNTERCLOCKWISE));
		else
			return polygon.is_simple();
	}

	static bool isSimplePolygon(const std::vector<Complex>& polygon, bool bCheckOrientation = true);
	static bool isSimplePolygonInexact(const std::vector<Complex>& polygon, bool bCheckOrientation = true);
	static bool triangulatePolygonSafe(const std::vector<std::complex<double> >& polygonPoints, std::vector<std::complex<double> >& meshVertices, std::vector<unsigned int>& triangleIndices, std::vector<std::pair<int, std::complex<double> > >& boundaryVertices, double maxTriangleArea, bool subsampleBoundaryEdges);
};

typedef CompilationAccelerator CA;
