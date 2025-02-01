#include <iostream>
#include <fstream>
#include <boost/format.hpp>

#include "triangulation.h"
#include "cassert"

extern "C"
{
	#include "triangle.h"
}

// ZL this is the version I'm working with; ignore the copy-pasted version afterwards.
bool triangulatePolygon(const std::vector<std::complex<double> >& polygonPoints, std::vector<std::complex<double> >& meshVertices, std::vector<unsigned int>& triangleIndices, std::vector<std::pair<int, std::complex<double> > >& boundaryVertices, double maxTriangleArea, bool subsampleBoundaryEdges, const char *triangle_params, std::vector<int> *additional_segments)
{
	int n = polygonPoints.size();

	if (n < 3)
	{
		return false;
	}

	struct triangulateio in, out;
	memset(&in, 0, sizeof(triangulateio));
	memset(&out, 0, sizeof(triangulateio));
	char refine_switchesAC[64] = { 0 };

	int n_additional_segments = additional_segments ? additional_segments->size() / 2 : 0;
	int nsegments = n_additional_segments + n;
	//printf("%d points and segments, %d additional (inner) segments\n", n, n_additional_segments);

	in.numberofpoints = n;
	in.numberofsegments = nsegments;
	in.numberofpointattributes = 1;
	in.pointlist = (double*)malloc((size_t)(2 * n * sizeof(double)));
	in.pointmarkerlist = (int*)malloc((size_t)(n * sizeof(int)));;
	in.segmentlist = (int*)malloc((size_t)(2 * nsegments * sizeof(int)));
	in.segmentmarkerlist = (int*)malloc((size_t)(nsegments * sizeof(int)));

	memcpy(in.pointlist, &(polygonPoints.front()), 2 * n * sizeof(double));

	for (int i = 0; i < n; i++)
	{
		in.pointmarkerlist[i] = 1;

		int currentIndex = i;
		int nextIndex = (i + 1 >= n) ? 0 : i + 1;

		in.segmentlist[2 * i] = currentIndex;
		in.segmentlist[2 * i + 1] = nextIndex;
		in.segmentmarkerlist[i] = 1;
	}

	// additional segments
	if (additional_segments && !additional_segments->empty()) {
		memcpy(in.segmentlist + 2 * n, &additional_segments->at(0), additional_segments->size() * sizeof(int));
		for (int i = 0; i < n_additional_segments; ++i)
			in.segmentmarkerlist[n + i] = 0;
	}

	// debug - save .poly
	// If command line triangle works fine, then it might be a precision issue. Round the coordinates to 1e-6 precision for a more accurate comparison.
	if (0) {
		std::ofstream fl("c:\\8\\disc.poly");
		fl << boost::format("%d 2 1 0\n") % n;
		for (int i = 0; i < n; i++)
			fl << boost::format("%d %f %f %d\n") % i % polygonPoints[i].real() % polygonPoints[i].imag() % in.pointmarkerlist[i];

		fl << boost::format("\n%d 1\n") % nsegments;
		for (int i = 0; i < nsegments; i++)
			fl << boost::format("%d %d %d %d\n") % i % in.segmentlist[2 * i] % in.segmentlist[2 * i + 1] % in.segmentmarkerlist[i];

		fl << "\n0\n";
	}

	in.numberoftriangles = 0;
	in.numberoftriangleattributes = 0;
	in.numberofpointattributes = 0;
	in.numberofcorners = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.regionlist = NULL;

	int minAngle = 33; // 33

	if (triangle_params)
		strncpy(refine_switchesAC, triangle_params, 64);
	else {
		if (subsampleBoundaryEdges)
		{
			sprintf_s(refine_switchesAC, 64, "Qpzq%da%f", minAngle, maxTriangleArea);
		}
		else
		{
			sprintf_s(refine_switchesAC, 64, "YQpzq%da%f", minAngle, maxTriangleArea);
		}
	}
	//refine_switchesAC[0] = 0; // test with default params

	//std::cout << "Calling triangle with " << refine_switchesAC << std::endl;

	triangulate(refine_switchesAC, &in, &out, NULL);

	if (out.numberoftriangles <= 0 || out.numberofsegments <= 0 || out.numberofpoints <= 0)
	{
		std::cerr << "ERROR in triangulatePolygon: triangulation failed" << std::endl;
		return false;
	}
	else {
		//std::cout << "triangle finished" << std::endl;
	}

	int p = out.numberofpoints;
	meshVertices.resize(p);
	boundaryVertices.clear();
	for (int i = 0; i < p; i++)
	{
		std::complex<double> position(out.pointlist[2 * i], out.pointlist[2 * i + 1]);
		meshVertices[i] = position;
		//std::cout << i << ": " << position << std::endl;
		if (out.pointmarkerlist[i] == 1)
		{
			boundaryVertices.push_back(std::pair<int, std::complex<double> >(i, position));
		}
	}

	int t = out.numberoftriangles;
	triangleIndices.resize(3 * t);
	for (int i = 0; i < t; i++)
	{
		assert(out.trianglelist[3 * i] >= 0);
		assert(out.trianglelist[3 * i + 1] >= 0);
		assert(out.trianglelist[3 * i + 2] >= 0);
		triangleIndices[3 * i] = out.trianglelist[3 * i];
		triangleIndices[3 * i + 1] = out.trianglelist[3 * i + 1];
		triangleIndices[3 * i + 2] = out.trianglelist[3 * i + 2];
		//printf("t%d: %d %d %d\n", i, triangleIndices[3*i], triangleIndices[3*i+1], triangleIndices[3*i+2]);
	}

	freeTriangleStructure(in);
	freeTriangleStructure(out);

	return true;
}
bool triangulatePolygonWithHoles(const std::vector<std::complex<double> >& polygonPoints, const std::vector<std::complex<double> >& holesCoordinatesList, std::vector<std::complex<double> >& meshVertices, std::vector<unsigned int>& triangleIndices, std::vector<std::pair<int, std::complex<double> > >& boundaryVertices, double maxTriangleArea, bool subsampleBoundaryEdges, const char *triangle_params, std::vector<int> *additional_segments)
	{
	int n = polygonPoints.size();
	int n1 = holesCoordinatesList.size();

	if(n < 3)
	{
		return false;
	}

	struct triangulateio in, out;
	memset(&in, 0, sizeof(triangulateio));
	memset(&out, 0, sizeof(triangulateio));
	char refine_switchesAC[64] = {0};

    int n_additional_segments = additional_segments ? additional_segments->size() / 2 : 0;
    int nsegments = n_additional_segments + n;
//printf("%d points and segments, %d additional (inner) segments\n", n, n_additional_segments);

	in.numberofpoints = n;
	in.numberofholes = n1;
	in.numberofsegments = nsegments;
	in.numberofpointattributes = 1;
	in.pointlist = (double*)malloc((size_t)(2*n*sizeof(double)));
	in.holelist = (double*)malloc((size_t)(2*n1*sizeof(double)));
	in.pointmarkerlist = (int*)malloc((size_t)(n*sizeof(int)));;
	in.segmentlist = (int*)malloc((size_t)(2*nsegments*sizeof(int)));
	in.segmentmarkerlist = (int*)malloc((size_t)(nsegments*sizeof(int)));

	memcpy(in.pointlist, &(polygonPoints.front()), 2*n*sizeof(double));
	memcpy(in.holelist, &(holesCoordinatesList.front()), 2*n1*sizeof(double));

	for(int i = 0; i < n; i++)
	{
		in.pointmarkerlist[i] = 1;

		int currentIndex = i;
		int nextIndex = (i+1 >= n) ? 0 : i+1;

		in.segmentlist[2*i] = currentIndex;
		in.segmentlist[2*i+1] = nextIndex;
		in.segmentmarkerlist[i] = 1;
	}       

    // additional segments
    if ( additional_segments && !additional_segments->empty() ) {
        memcpy(in.segmentlist + 2*n, &additional_segments->at(0), additional_segments->size()*sizeof(int));
        for ( int i = 0; i < n_additional_segments; ++i ) 
            in.segmentmarkerlist[n + i] = 0;
	}

    // debug - save .poly
    // If command line triangle works fine, then it might be a precision issue. Round the coordinates to 1e-6 precision for a more accurate comparison.
    if ( 0 ) {
		std::ofstream fl("c:\\8\\disc.poly");
        fl << boost::format("%d 2 1 0\n") % n;
        for ( int i = 0; i < n; i++ )
            fl << boost::format("%d %f %f %d\n") % i % polygonPoints[i].real() % polygonPoints[i].imag() % in.pointmarkerlist[i];

		fl << boost::format("\n%d 1\n") % nsegments;
		for ( int i = 0; i < nsegments; i++ )
			fl << boost::format("%d %d %d %d\n") % i % in.segmentlist[2*i] % in.segmentlist[2*i+1] % in.segmentmarkerlist[i];

        fl << "\n0\n";
	}

	in.numberoftriangles = 0;
	in.numberoftriangleattributes = 0;
	in.numberofpointattributes = 0;
	in.numberofcorners = 0;
	in.numberofregions = 0;
	in.regionlist = NULL;
	int minAngle = 0.0000001; // 33 was 0.1
	
	if ( triangle_params )
		strncpy(refine_switchesAC, triangle_params, 64);
	else {
		if(subsampleBoundaryEdges)
		{
			sprintf_s(refine_switchesAC, 64, "Qpzq%da%f", minAngle, maxTriangleArea);		
		}
		else
		{
			sprintf_s(refine_switchesAC, 64, "YQpzq%da%f", minAngle, maxTriangleArea);
		}
	}
	//refine_switchesAC[0] = 0; // test with default params

	//std::cout << "Calling triangle with " << refine_switchesAC << std::endl;

	triangulate(refine_switchesAC, &in, &out, NULL);

	if(out.numberoftriangles <= 0 || out.numberofsegments <= 0 || out.numberofpoints <= 0)
	{
		std::cerr << "ERROR in triangulatePolygon: triangulation failed" << std::endl;
		return false;
	} else {
		//std::cout << "triangle finished" << std::endl;
	}

	int p = out.numberofpoints;
	meshVertices.resize(p);
	boundaryVertices.clear();
	for(int i = 0; i < p; i++)
	{
		std::complex<double> position(out.pointlist[2*i], out.pointlist[2*i + 1]);
		meshVertices[i] = position;
//std::cout << i << ": " << position << std::endl;
		if(out.pointmarkerlist[i] == 1)
		{
			boundaryVertices.push_back(std::pair<int, std::complex<double> >(i, position));
		}
	}

	int t = out.numberoftriangles;
	triangleIndices.resize(3*t);
	for(int i = 0; i < t; i++)
	{
		assert(out.trianglelist[3*i] >= 0);
		assert(out.trianglelist[3*i+1] >= 0);
		assert(out.trianglelist[3*i+2] >= 0);
		triangleIndices[3*i] = out.trianglelist[3*i];
		triangleIndices[3*i+1] = out.trianglelist[3*i+1];
		triangleIndices[3*i+2] = out.trianglelist[3*i+2];
//printf("t%d: %d %d %d\n", i, triangleIndices[3*i], triangleIndices[3*i+1], triangleIndices[3*i+2]);
	}

	freeTriangleStructure(in);
	freeTriangleStructure(out);

	return true;
}


bool triangulatePolygonWithoutAddingVertices(const std::vector<std::complex<double> >& polygonPoints, std::vector<unsigned int>& triangleIndices)
{
	int n = polygonPoints.size();

	if(n < 3)
	{
		return false;
	}
	for(int i = 0; i < n; i++)
	{
		int nextIndex = i < n-1 ? i+1 : 0;
		if(polygonPoints[i] == polygonPoints[nextIndex])
		{
			std::cerr << "Cannot triangulate degenerated polygon!\n";
			return false;
		}
	}
	struct triangulateio in, out;
	memset(&in, 0, sizeof(triangulateio));
	memset(&out, 0, sizeof(triangulateio));

	in.numberofpoints = n;
	in.numberofsegments = n;
	in.numberofpointattributes = 1;
	in.pointlist = (double*)malloc((size_t)(2*n*sizeof(double)));
	in.segmentlist = (int*)malloc((size_t)(2*n*sizeof(int)));
	in.segmentmarkerlist = (int*)malloc((size_t)(n*sizeof(int)));

	memcpy(in.pointlist, &(polygonPoints.front()), 2*n*sizeof(double));

	for(int i = 0; i < n; i++)
	{
		int currentIndex = i;
		int nextIndex = (i+1 >= n) ? 0 : i+1;

		in.segmentlist[2*i] = currentIndex;
		in.segmentlist[2*i+1] = nextIndex;
		in.segmentmarkerlist[i] = 1;
	}       

	in.pointmarkerlist = NULL;
	in.numberoftriangles = 0;
	in.numberoftriangleattributes = 0;
	in.numberofpointattributes = 0;
	in.numberofcorners = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.regionlist = NULL;

	char refine_switchesAC[] = "QpYYz";
	triangulate(refine_switchesAC, &in, &out, NULL);

	if(out.numberoftriangles <= 0 || out.numberofsegments <= 0 || out.numberofpoints <= 0)
	{
		std::cerr << "ERROR in triangulatePolygon: triangulation failed" << std::endl;
		return false;
	}
	int t = out.numberoftriangles;
	triangleIndices.resize(3*t);

	for(int i = 0; i < t; i++)
	{
		assert(out.trianglelist[3*i] >= 0);
		assert(out.trianglelist[3*i+1] >= 0);
		assert(out.trianglelist[3*i+2] >= 0);
		triangleIndices[3*i] = out.trianglelist[3*i];
		triangleIndices[3*i+1] = out.trianglelist[3*i+1];
		triangleIndices[3*i+2] = out.trianglelist[3*i+2];
	}

	freeTriangleStructure(in);
	freeTriangleStructure(out);

	return true;
}
bool triangulatePolygonWithHolesWithoutAddingVertices(const std::vector<std::complex<double> >& polygonPoints, const std::vector<std::complex<double> >& holesCoordinatesList, int numInternalPoints, std::vector<unsigned int>& triangleIndices)
{
	int n = polygonPoints.size();
	int n1 = holesCoordinatesList.size();
	int numOuterPoints = n - numInternalPoints;

	if (n < 3)
	{
		return false;
	}

	struct triangulateio in, out;
	memset(&in, 0, sizeof(triangulateio));
	memset(&out, 0, sizeof(triangulateio));
	char refine_switchesAC[] = "QpYYz";

	/*int n_additional_segments = additional_segments ? additional_segments->size() / 2 : 0;
	int nsegments = n;*/
	//printf("%d points and segments, %d additional (inner) segments\n", n, n_additional_segments);

	in.numberofpoints = n;
	in.numberofholes = n1;
	in.numberofsegments = n;
	in.numberofpointattributes = 1;
	in.pointlist = (double*)malloc((size_t)(2 * n * sizeof(double)));
	in.holelist = (double*)malloc((size_t)(2 * n1 * sizeof(double)));
	in.pointmarkerlist = (int*)malloc((size_t)(n * sizeof(int)));;
	in.segmentlist = (int*)malloc((size_t)(2 * n * sizeof(int)));
	in.segmentmarkerlist = (int*)malloc((size_t)(n * sizeof(int)));

	memcpy(in.pointlist, &(polygonPoints.front()), 2 * n * sizeof(double));
	memcpy(in.holelist, &(holesCoordinatesList.front()), 2 * n1 * sizeof(double));

	for (int i = 0; i < n; i++)
	{
		in.pointmarkerlist[i] = 1;

		int currentIndex = i;
		int nextIndex = 0;
		//int nextIndex = (i + 1 >= n) ? 0 : i + 1;
		if (i < numOuterPoints)
		{
			nextIndex = (i + 1 >= numOuterPoints) ? 0 : i + 1;
		}
		else
		{
			nextIndex = (i + 1 >= n) ? numOuterPoints : i + 1;
		}


		in.segmentlist[2 * i] = currentIndex;
		in.segmentlist[2 * i + 1] = nextIndex;
		in.segmentmarkerlist[i] = 1;
	}

	// additional segments
	/*if (additional_segments && !additional_segments->empty()) {
		memcpy(in.segmentlist + 2 * n, &additional_segments->at(0), additional_segments->size() * sizeof(int));
		for (int i = 0; i < n_additional_segments; ++i)
			in.segmentmarkerlist[n + i] = 0;
	}*/

	in.numberoftriangles = 0;
	in.numberoftriangleattributes = 0;
	in.numberofpointattributes = 0;
	in.numberofcorners = 0;
	in.numberofregions = 0;
	in.regionlist = NULL;
	int minAngle = 33; // 33

	//refine_switchesAC[0] = 0; // test with default params

	//std::cout << "Calling triangle with " << refine_switchesAC << std::endl;

	triangulate(refine_switchesAC, &in, &out, NULL);

	if (out.numberoftriangles <= 0 || out.numberofsegments <= 0 || out.numberofpoints <= 0)
	{
		std::cerr << "ERROR in triangulatePolygon: triangulation failed" << std::endl;
		return false;
	}
	else {
		//std::cout << "triangle finished" << std::endl;
	}

	int p = out.numberofpoints;
	int t = out.numberoftriangles;
	triangleIndices.resize(3 * t);
	for (int i = 0; i < t; i++)
	{
		assert(out.trianglelist[3 * i] >= 0);
		assert(out.trianglelist[3 * i + 1] >= 0);
		assert(out.trianglelist[3 * i + 2] >= 0);
		triangleIndices[3 * i] = out.trianglelist[3 * i];
		triangleIndices[3 * i + 1] = out.trianglelist[3 * i + 1];
		triangleIndices[3 * i + 2] = out.trianglelist[3 * i + 2];
		//printf("t%d: %d %d %d\n", i, triangleIndices[3*i], triangleIndices[3*i+1], triangleIndices[3*i+2]);
	}

	freeTriangleStructure(in);
	freeTriangleStructure(out);

	return true;

}


void freeTriangleStructure(struct triangulateio& triangleStruct)
{
	if(triangleStruct.pointlist)
	{
		free(triangleStruct.pointlist);
		triangleStruct.pointlist = NULL;
	}
	if(triangleStruct.pointattributelist)
	{
		free(triangleStruct.pointattributelist);
		triangleStruct.pointattributelist = NULL;
	}
	if(triangleStruct.pointmarkerlist)
	{
		free(triangleStruct.pointmarkerlist);
		triangleStruct.pointmarkerlist = NULL;
	}
	if(triangleStruct.regionlist)
	{
		free(triangleStruct.regionlist);
		triangleStruct.regionlist = NULL;
	}
	if(triangleStruct.pointlist)
	{
		free(triangleStruct.pointlist);
		triangleStruct.pointlist = NULL;
	}

	if(triangleStruct.pointattributelist)
	{
		free(triangleStruct.pointattributelist);
		triangleStruct.pointattributelist = NULL;
	}
	if(triangleStruct.trianglelist)
	{
		free(triangleStruct.trianglelist);
		triangleStruct.trianglelist = NULL;
	}
	if(triangleStruct.triangleattributelist)
	{
		free(triangleStruct.triangleattributelist);
		triangleStruct.triangleattributelist = NULL;
	}
	if(triangleStruct.trianglearealist)
	{
		free(triangleStruct.trianglearealist);
		triangleStruct.trianglearealist = NULL;
	}
	if(triangleStruct.segmentlist)
	{
		free(triangleStruct.segmentlist);
		triangleStruct.segmentlist = NULL;
	}
	if(triangleStruct.segmentmarkerlist)
	{
		free(triangleStruct.segmentmarkerlist);
		triangleStruct.segmentmarkerlist = NULL;
	}
	if(triangleStruct.edgelist)
	{
		free(triangleStruct.edgelist);
		triangleStruct.edgelist = NULL;
	}
	if(triangleStruct.edgemarkerlist)
	{
		free(triangleStruct.edgemarkerlist);
		triangleStruct.edgemarkerlist = NULL;
	}
}



void Exactinit()
{
	exactinit();
}



double CounterClockWise(const double pa[2], const double pb[2], const double pc[2])
{
	return orient2D((double*)pa, (double*)pb, (double*)pc);
}



double CounterClockWise(const std::complex<double>& pa, const std::complex<double>& pb, const std::complex<double>& pc)
{
	return orient2D((double*)(&pa), (double*)(&pb), (double*)(&pc));
}
