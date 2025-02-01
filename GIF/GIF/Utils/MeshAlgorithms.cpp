#include "stdafx.h"

#include "MeshAlgorithms.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "CGAL/CGAL_Macros.h"
#include "CGAL/CGAL_UTIL.h"
#include "Utils/Utilities.h"
#include "Utils/triangulation.h"
void MeshAlgorithms::buildTriangleMesh(const std::vector<Point_3>& vertices, const std::vector<unsigned int>& faceIndices, Mesh& mesh)
{
	//we assume here that the indices are not negative!
	MeshBuilder<HDS, Kernel> meshBuilder(&vertices, (std::vector<int>*)(&faceIndices), NULL);
	mesh.clear();
	mesh.delegate(meshBuilder);
	mesh.updateAllGlobalIndices();
}


int MeshAlgorithms::triangulatePolygonWithHolesAndCreatedCgalMesh(std::vector<Complex>& polygon, const std::vector<std::complex<double> >& holesCoordinatesList, int numWantedTriangles, bool subsampleBoundary, Mesh& mesh, const char *triangle_params, std::vector<int> *additional_segments)
{
	double polygonTotalArea = computePolygonArea(polygon);
	double averageTriangleArea = polygonTotalArea / numWantedTriangles;
	//we cannot feed the "triangle Shewchuk code" with average triangle area. instead, we can feed it with maximum triangle area.
	//we use this conversion between average and maximum area based on experiments
	double maximumTriangleArea = 1.7*averageTriangleArea;

	std::vector<Complex> meshVertices;
	std::vector<std::pair<int, Complex> > boundaryVertices;
	std::vector<unsigned int> triangleIndices;

	bool result = triangulatePolygonWithHoles(polygon, holesCoordinatesList, meshVertices, triangleIndices, boundaryVertices, maximumTriangleArea, subsampleBoundary, triangle_params, additional_segments);

	if (!result)
	{
		return 0;
	}

	int numV = meshVertices.size();
	std::vector<Point_3> vertices(numV);

	for (int i = 0; i < numV; i++)
	{
		vertices[i] = Point_3(CGUT::complexToPoint3(meshVertices[i]));
	}

	MeshAlgorithms::buildTriangleMesh(vertices, triangleIndices, mesh);

	assert(mesh.is_pure_triangle());

	int numTriangles = mesh.size_of_facets();

	return numTriangles;
}

int MeshAlgorithms::triangulatePolygonWithHolesWithoutAddingVerticesAndCreatedCgalMesh(std::vector<Complex>& polygon, const std::vector<std::complex<double> >& holesCoordinatesList, int numInternalPoints, Mesh& mesh)
{
	/*double polygonTotalArea = computePolygonArea(polygon);
	double averageTriangleArea = polygonTotalArea / numWantedTriangles;*/
	//we cannot feed the "triangle Shewchuk code" with average triangle area. instead, we can feed it with maximum triangle area.
	//we use this conversion between average and maximum area based on experiments
	//double maximumTriangleArea = 1.7*averageTriangleArea;

	/*std::vector<Complex> meshVertices;
	std::vector<std::pair<int, Complex> > boundaryVertices;*/
	std::vector<unsigned int> triangleIndices;

	bool result = triangulatePolygonWithHolesWithoutAddingVertices(polygon, holesCoordinatesList, numInternalPoints, triangleIndices);

	if (!result)
	{
		return 0;
	}

	GMMDenseColMatrix scaffoldBorderPoints(polygon.size(), 2);
	for (int i = 0; i < polygon.size(); i++)
	{
		scaffoldBorderPoints(i, 0) = real(polygon[i]);
		scaffoldBorderPoints(i, 1) = imag(polygon[i]);
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("borderPoints1", scaffoldBorderPoints);

	int numV = polygon.size();
	std::vector<Point_3> vertices(numV);

	for (int i = 0; i < numV; i++)
	{
		vertices[i] = Point_3(CGUT::complexToPoint3(polygon[i]));
	}

	MeshAlgorithms::buildTriangleMesh(vertices, triangleIndices, mesh);

	assert(mesh.is_pure_triangle());

	int numTriangles = mesh.size_of_facets();

	return numTriangles;
}
int MeshAlgorithms::buildScaffoldWithFreeOuterBoundary(std::vector<std::complex<double>> internalBorderOfScaffold, bool addVertices, Mesh& scaffold, int sizeOfOuterBoundary, double radiusesRatio, bool isNarrow)
{
	int numOfTriangles = 350;
	if (addVertices)
	{
		int i = 0, r = internalBorderOfScaffold.size();
		double x1 = 0.0, y1 = 0.0;
		/*std::vector<std::complex<double>> internalBorderOfScaffold(r);
		for (i = 0; i < r; i++)
		{
			internalBorderOfScaffold[i] = Complex(internalBounary(i, 0), internalBounary(i, 1));
		}*/
		i = 0;
		for (i = 0; i < r; i++)
		{
			x1 += internalBorderOfScaffold[i].real() / r;
			y1 += internalBorderOfScaffold[i].imag() / r;
		}//(x1, y1) is the average value of the u and v coordinates of all the border vertices from the cgalMesh's uv set

		double radiusOfBigCircle = 0.0;
		for (i = 0; i < r; i++)
		{
			double curDis = sqrt((internalBorderOfScaffold[i].real() - x1) * (internalBorderOfScaffold[i].real() - x1) + (internalBorderOfScaffold[i].imag() - y1) * (internalBorderOfScaffold[i].imag() - y1));
			if (curDis > radiusOfBigCircle)
				radiusOfBigCircle = curDis;
		}

		double r1 = radiusOfBigCircle;
		std::vector<int> additional_segments[1];
		double angle = 0.0;
		radiusOfBigCircle *= radiusesRatio;//the radius of the inner circle equals to the distance between (x, y) and the farthest border point from it, multiplied by 2
		if (sizeOfOuterBoundary == 1)
			r *= 4;
		else
			r = sizeOfOuterBoundary;
		std::vector<std::complex<double>> allPolygonPoints(r + internalBorderOfScaffold.size());//all the border points of the annulus
		angle = 0.0;
		for (i = 0; i < r; i++)//enter the values of the external border points of the annulus
		{
			allPolygonPoints[i] = Complex(radiusOfBigCircle * cos(angle), radiusOfBigCircle * sin(angle));
			angle += 2.0 * M_PI / (double)r;
			if (i != r - 1)
			{
				additional_segments[0].push_back(i);
				additional_segments[0].push_back(i + 1);
			}
			else
			{
				additional_segments[0].push_back(r - 1);
				additional_segments[0].push_back(0);
			}
		}
		for (i = 0; i < internalBorderOfScaffold.size(); i++)//enter the values of the internal border points of the annulus
		{
			allPolygonPoints[r + i] = internalBorderOfScaffold[i];
			if (i != internalBorderOfScaffold.size() - 1)
			{
				additional_segments[0].push_back(r + i);
				additional_segments[0].push_back(r + i + 1);
			}
			else
			{
				additional_segments[0].push_back(r + internalBorderOfScaffold.size() - 1);
				additional_segments[0].push_back(r);
			}
		}
		/*cout << additional_segments[0] << endl;
		cout << allPolygonPoints << endl;*/
		std::vector<std::complex<double> > holesCoordinatesList(1);
		holesCoordinatesList[0] = Complex(x1, y1);
		bool result = MeshAlgorithms::triangulatePolygonWithHolesAndCreatedCgalMesh(allPolygonPoints, holesCoordinatesList, numOfTriangles, false, scaffold, nullptr, additional_segments);
	}
	else
	{
		int i = 0, r = internalBorderOfScaffold.size();
		double x1 = 0.0, y1 = 0.0;
		/*std::vector<std::complex<double>> internalBorderOfScaffold(r);
		for (i = 0; i < r; i++)
		{
			internalBorderOfScaffold[i] = Complex(internalBounary(i, 0), internalBounary(i, 1));
		}*/
		i = 0;
		for (i = 0; i < r; i++)
		{
			x1 += internalBorderOfScaffold[i].real() / r;
			y1 += internalBorderOfScaffold[i].imag() / r;
		}//(x1, y1) is the average value of the u and v coordinates of all the border vertices from the cgalMesh's uv set

		double radiusOfBigCircle = 0.0;
		for (i = 0; i < r; i++)
		{
			double curDis = sqrt((internalBorderOfScaffold[i].real() - x1) * (internalBorderOfScaffold[i].real() - x1) + (internalBorderOfScaffold[i].imag() - y1) * (internalBorderOfScaffold[i].imag() - y1));
			if (curDis > radiusOfBigCircle)
				radiusOfBigCircle = curDis;
		}

		double r1 = radiusOfBigCircle;
		std::vector<int> additional_segments[1];
		double angle = 0.0;
		radiusOfBigCircle *= radiusesRatio;//the radius of the inner circle equals to the distance between (x, y) and the farthest border point from it, multiplied by 2
		if (sizeOfOuterBoundary == 1)
			r *= 4;
		else
			r = sizeOfOuterBoundary;
		std::vector<std::complex<double>> allPolygonPoints(r + internalBorderOfScaffold.size());//all the border points of the annulus
		angle = 0.0;
		for (i = 0; i < r; i++)//enter the values of the external border points of the annulus
		{
			allPolygonPoints[i] = Complex(radiusOfBigCircle * cos(angle), radiusOfBigCircle * sin(angle));
			angle += 2.0 * M_PI / (double)r;
		}
		for (i = 0; i < internalBorderOfScaffold.size(); i++)//enter the values of the internal border points of the annulus
		{
			allPolygonPoints[r + i] = internalBorderOfScaffold[i];
		}
		/*cout << additional_segments[0] << endl;
		cout << allPolygonPoints << endl;*/
		std::vector<std::complex<double> > holesCoordinatesList(1);
		holesCoordinatesList[0] = Complex(x1, y1);
		int num1 = internalBorderOfScaffold.size();
		if (!isNarrow)
			bool result = MeshAlgorithms::triangulatePolygonWithHolesWithoutAddingVerticesAndCreatedCgalMesh(allPolygonPoints, holesCoordinatesList, num1, scaffold);
		else
		{
			GMMDenseColMatrix scaffoldInternalBorderPoints(internalBorderOfScaffold.size(), 2);
			GMMDenseColMatrix scaffoldExternalBorderPoints(r, 2);
			for (int i = 0; i < internalBorderOfScaffold.size(); i++)
			{
				scaffoldInternalBorderPoints(i, 0) = real(internalBorderOfScaffold[i]);
				scaffoldInternalBorderPoints(i, 1) = imag(internalBorderOfScaffold[i]);
			}
			for (int i = 0; i < r; i++)
			{
				scaffoldExternalBorderPoints(i, 0) = real(allPolygonPoints[i]);
				scaffoldExternalBorderPoints(i, 1) = imag(allPolygonPoints[i]);
			}
			MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldInternalBorderPoints1", scaffoldInternalBorderPoints);
			MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldExternalBorderPoints", scaffoldExternalBorderPoints);

		}
	}
	return 0;
}

int MeshAlgorithms::buildScaffoldWithFixedOuterBoundary(std::vector<std::complex<double>> allPolygonPoints, bool addVertices, Mesh& scaffold, int numOfInternalVertices)
{
	int numOfTriangles = 350;
	if (addVertices)
	{
		//std::vector<std::complex<double>> internalPolygonPoints1 = internalPolygonPoints[0]; //all the scaffold vertices
		int i = 0, r = allPolygonPoints.size() - numOfInternalVertices;
		std::vector<int> additional_segments[1];
		for (i = 0; i < r; i++)//enter the values of the external border points of the annulus
		{
			if (i != r - 1)
			{
				additional_segments[0].push_back(i);
				additional_segments[0].push_back(i + 1);
			}
			else
			{
				additional_segments[0].push_back(r - 1);
				additional_segments[0].push_back(0);
			}
		}
		for (i = 0; i < numOfInternalVertices; i++)//enter the values of the internal border points of the annulus
		{
			if (i != numOfInternalVertices - 1)
			{
				additional_segments[0].push_back(r + i);
				additional_segments[0].push_back(r + i + 1);
			}
			else
			{
				additional_segments[0].push_back(r + numOfInternalVertices - 1);
				additional_segments[0].push_back(r);
			}
		}
		/*cout << additional_segments[0] << endl;
		cout << internalPolygonPoints1 << endl;*/
		std::vector<std::complex<double> > holesCoordinatesList(1);
		/*Complex next = internalBorderOfScaffold[2];
		Complex prev = internalBorderOfScaffold[0];
		Complex cur = internalBorderOfScaffold[1];*/

		/*Complex firstRatio = (next - prev) / (cur - prev);
		firstRatio = sqrt(firstRatio);
		Complex secondRatio = (cur - next) / (prev - next);
		secondRatio = sqrt(secondRatio);
		Complex F = (secondRatio / firstRatio) * Complex(0, -1) * (next - prev);
		holesCoordinatesList[0] = cur + 0.001 * F / norm(F);*/
		GMMDenseColMatrix holeCoordinates(2, 1);
		MatlabGMMDataExchange::GetEngineDenseMatrix("holeCoor", holeCoordinates);
		holesCoordinatesList[0] = Complex(holeCoordinates(0, 0), holeCoordinates(0, 1));
		/*GMMDenseColMatrix scaffoldBorderPoints(internalPolygonPoints1.size(), 2);
		for (int i = 0; i < internalPolygonPoints1.size(); i++)
		{
			scaffoldBorderPoints(i, 0) = real(internalPolygonPoints1[i]);
			scaffoldBorderPoints(i, 1) = imag(internalPolygonPoints1[i]);
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldBorderPoints1", scaffoldBorderPoints);*/
		//holesCoordinatesList[0] = Complex(-0.6, 0.0);
		//int num1 = internalBorderOfScaffold.size();
		bool result = MeshAlgorithms::triangulatePolygonWithHolesAndCreatedCgalMesh(allPolygonPoints, holesCoordinatesList, numOfTriangles, false, scaffold, nullptr, additional_segments);
	}
	else
	{

		//std::vector<std::complex<double>> internalPolygonPoints1 = internalPolygonPoints[0]; //all the scaffold vertices
		/*GMMDenseColMatrix holeCoordinates(2, 1);
		MatlabGMMDataExchange::GetEngineDenseMatrix("holeCoor", holeCoordinates);*/
		//holesCoordinatesList[0] = Complex(holeCoordinates(0, 0), holeCoordinates(0, 1));
		/*GMMDenseColMatrix scaffoldBorderPoints(internalPolygonPoints1.size(), 2);
		for (int i = 0; i < internalPolygonPoints1.size(); i++)
		{
			scaffoldBorderPoints(i, 0) = real(internalPolygonPoints1[i]);
			scaffoldBorderPoints(i, 1) = imag(internalPolygonPoints1[i]);
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldBorderPoints1", scaffoldBorderPoints);*/
		//holesCoordinatesList[0] = Complex(-0.6, 0.0);
		std::vector<std::complex<double> > holesCoordinatesList(1);
		GMMDenseColMatrix holeCoordinates(2, 1);
		MatlabGMMDataExchange::GetEngineDenseMatrix("holeCoor", holeCoordinates);
		holesCoordinatesList[0] = Complex(holeCoordinates(0, 0), holeCoordinates(0, 1));
		//int num1 = internalBorderOfScaffold.size();
		bool result = MeshAlgorithms::triangulatePolygonWithHolesWithoutAddingVerticesAndCreatedCgalMesh(allPolygonPoints, holesCoordinatesList, numOfInternalVertices, scaffold);
	}
	return 0;
}