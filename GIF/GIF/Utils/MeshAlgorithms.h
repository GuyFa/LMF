
#pragma once

#include "CGAL/CGAL_Mesh.h"
#include "CGAL/CGAL_Macros.h"
#include "Utils/GMM_Macros.h"
#include "utils/STL_Macros.h"

typedef complex<double> Complex;

namespace MeshAlgorithms
{
	void buildTriangleMesh(const std::vector<Point_3>& vertices, const std::vector<unsigned int>& faceIndices, Mesh& mesh);
	int buildScaffoldWithFreeOuterBoundary(std::vector<std::complex<double>> internalBorderOfScaffold, bool addVertices, Mesh& scaffold, int sizeOfOuterBoundary = 1, double radiusesRatio = 2.0, bool isNorrow = false); //this function needs the internal vertices of the scaffold, and computes the outer boundary by itself (the outer boundary radius can be changed by setting the sizeOfOuterBoundary variable)
	int buildScaffoldWithFixedOuterBoundary(std::vector<std::complex<double>> allPolygonPoints, bool addVertices, Mesh& scaffold, int numOfInternalVertices); //this function needs all the vertices of the scaffold and the number of internal vertices, since all the the boundaries are already fixed (the variable allPolygonPoints contains all the vertices, where the outer boundary vertices are placed first)
	int triangulatePolygonWithHolesAndCreatedCgalMesh(std::vector<Complex>& polygon, const std::vector<std::complex<double> >& holesCoordinatesList, int numWantedTriangles, bool subsampleBoundary, Mesh& mesh, const char *triangle_params = nullptr, std::vector<int> *additional_segments = nullptr);
	int triangulatePolygonWithHolesWithoutAddingVerticesAndCreatedCgalMesh(std::vector<Complex>& polygon, const std::vector<std::complex<double> >& holesCoordinatesList, int numInternalPoints, Mesh& mesh);
}
