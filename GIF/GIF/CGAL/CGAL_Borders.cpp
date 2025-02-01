#include "Utils/stdafx.h"

#include "CGAL_Mesh.h"
#include "CGAL_Macros.h"
#include "CGAL_Borders.h"



CGAL_Borders::CGAL_Borders(Mesh& mesh) : mV(0), mF(0), mE(0), mC(0), mB(0)
{
	mV = mesh.size_of_vertices();
	mF = mesh.size_of_facets();
	int numHalfedges = mesh.size_of_halfedges();
	assert((numHalfedges % 2) == 0);
	mE = numHalfedges / 2;
	mC = computeNumberOfConnectedComponents(mesh);
	computeBorders(mesh);
}


void CGAL_Borders::printMeshTopologicalStatistics()
{
	cout << "#Vertices: " << mV << "\n";
	cout << "#Faces: " << mF << "\n";
	cout << "#Edges: " << mE << "\n";
	cout << "Genus: " << genus() << "\n";
	cout << "Euler Number: " << eulerNumber() << "\n";
	cout << "#Components: " << mC << "\n";
	cout << "#Border Curves: " << mB << "\n";
	if (mB > 0)
	{
		int numBorderVertices = 0;
		for (int i = 0; i < mB; i++)
		{
			cout << "#Vertices on border number " << i << ": ";
			int numVerticesOnCurrentBorder = mBorders[i].size();
			cout << numVerticesOnCurrentBorder << "\n";
			numBorderVertices += numVerticesOnCurrentBorder;
		}
		cout << "Total Number of Border Vertices: " << numBorderVertices << "\n";
	}
	cout << endl;
}


int CGAL_Borders::genus() const
{
	int g = 0;
	int numerator = 2 * mC + mE - mV - mF - mB;
	assert(numerator % 2 == 0);
	return numerator / 2;
}


//computes the Euler characteristic number Chi
int CGAL_Borders::eulerNumber() const
{
	int euler = mV + mF - mE;
	return euler;
}


int CGAL_Borders::numVertices() const
{
	return mV;
}


int CGAL_Borders::numFaces() const
{
	return mF;
}


int CGAL_Borders::numEdges() const
{
	return mE;
}

//the connectivity is based on the dual graph and assumes that the input is a 2-manifold mesh without singular vertices nor edges.
int CGAL_Borders::numConnectedComponents() const
{
	return mC;
}


int CGAL_Borders::numBorders() const
{
	return mB;
}


void CGAL_Borders::computeBorders(Mesh& mesh)
{
	//bool isBorderNormalize = mesh.normalized_border_is_valid(true);
	int k = mesh.size_of_border_halfedges(); //the amount of elements to reserve

	for (Halfedge_iterator he_border = mesh.border_halfedges_begin(); he_border != mesh.halfedges_end(); ++he_border)
	{
		makeBorderHalfedge(he_border);
		he_border->userIndex() = 0; //mark as unvisited
	}

	mB = 0;
	mBorders.reserve(20);

	//The range [border_halfedges_begin(), halfedges_end()) denotes all border edges.
	//The halfedges are not necessarily grouped into different borders and inside each border, they are not necessarily ordered.
	Halfedge_iterator he_border = mesh.border_halfedges_begin();

	for (Halfedge_iterator he_border = mesh.border_halfedges_begin(); he_border != mesh.halfedges_end(); ++he_border)
	{
		makeBorderHalfedge(he_border);

		if (he_border->userIndex() == 0)
		{
			assert(k >= 2); //I guess the smallest boundary has at least 2 edges (it needs to be a hole for that).
			std::vector<Vertex_handle> verticesOfCurrentBorder;
			mBorders.push_back(verticesOfCurrentBorder);
			mBorders.back().reserve(k);
			int numVerticesInCurrentBorder = 0;
			do
			{
				//cout << he_border->vertex()->index() << endl; //debugging
				mBorders[mB].push_back(he_border->vertex());
				he_border->userIndex() = 1; //mark as visited
				he_border = he_border->next(); //go to the adjacent border half edge along the null face
				numVerticesInCurrentBorder++;
			} while (he_border->userIndex() == 0);
			k -= numVerticesInCurrentBorder;
			mB++;
		}
	}
}


//returns the number of connected components of the mesh by traversing the dual graph.
unsigned int CGAL_Borders::computeNumberOfConnectedComponents(Mesh& mesh)
{
	unsigned int numComponents = 0;

	for_each_facet(f, mesh)
	{
		f->user() = 0.0;
	}

	for_each_facet(f, mesh)
	{
		if (f->user() == 0.0)
		{
			numComponents++;
			traverseConnectedComponent(f, 0.0, 1.0);
		}
	}
	return numComponents;
}




//traversing the dual graph of the mesh starting from a given root face until all connected faces are reached.
//the function assumes that the unvisted faces are marked with user()==unVisitedTag and it marks every visited face by setting user() to visitedTag. 
void CGAL_Borders::traverseConnectedComponent(Facet_handle rootFacet, const double unVisitedTag, const double visitedTag)
{
	rootFacet->user() = visitedTag; //mark root as visited

	std::queue<Facet_handle> facetsQueue;
	facetsQueue.push(rootFacet);

	while (!facetsQueue.empty())
	{
		Facet_handle face = facetsQueue.front();
		facetsQueue.pop();

		Halfedge_around_facet_circulator h = face->facet_begin();
		Halfedge_around_facet_circulator h_end = h;

		CGAL_For_all(h, h_end)
		{
			Facet_handle neighborFace = h->opposite()->facet();
			if (neighborFace != NULL && neighborFace->user() == unVisitedTag)
			{
				neighborFace->user() = visitedTag;
				facetsQueue.push(neighborFace);
			}
		}
	}
}



//makes sure that the given halfedge points to the null face.
//the assumption is that the corresponding edge is on the border, meaning either the halfedge or its opposite are border.
void CGAL_Borders::makeBorderHalfedge(Halfedge_iterator& he_border)
{
	if (!he_border->is_border())
	{
		he_border = he_border->opposite(); //make sure the halfedge belongs to the null face. For some reason, this is not always the case
		assert(he_border->is_border()); //just to make sure it wasn't an interior halfedge.
	}
}


int CGAL_Borders::numVerticesInBorder(unsigned int index) const
{
	assert(index < mBorders.size());
	return mBorders[index].size();
}


bool CGAL_Borders::isTopologicalDisk() const
{
	if (mB == 1 && genus() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool CGAL_Borders::isMeshConnected() const
{
	return mC == 1;
}


Vertex_handle CGAL_Borders::vertex(unsigned int borderIndex, unsigned int vertexIndex) const
{
	assert(borderIndex < mBorders.size());
	assert(vertexIndex < mBorders[borderIndex].size());

	return mBorders[borderIndex][vertexIndex];
}


///////////////////////////////////////////////////////////
// computes the Geodesic curvature of the border curve //
///////////////////////////////////////////////////////////
// index - which border to use                         //
// useTargetMetric - whether to use the metric from the  //
//                   3D embedding or from the prescribed //
//                   edge lengths                        //
///////////////////////////////////////////////////////////
double CGAL_Borders::borderCurvature(unsigned int index, bool useTargetMetric) const
{
	assert(index < mBorders.size());

	double totalCurvature = 0.0;

	int numV = numVerticesInBorder(index);

	for (int i = 0; i < numV; i++)
	{
		Vertex_handle v = vertex(index, i);
		assert(v->is_border());

		double kappa = 0.0;

		if (useTargetMetric)
		{
			kappa = v->targetMetricGaussianCurvature();
		}
		else
		{
			kappa = v->gaussianCurvature();
		}

		totalCurvature += kappa;
	}
	return totalCurvature;
}

///////////////////////


//finds a list of border halfedges connecting "source" to "target"
//returns the uv length of the polyline defined by this path and 0.0 if there is no path connecting the two vertices
double CGAL_Borders::getBorderHalfEdgePath(const Vertex_handle& source, const Vertex_handle& target, std::vector<Halfedge_handle>& path)
{
	double length = 0.0;

	if (!source->is_border() || !target->is_border())
	{
		return 0.0;
	}

	Halfedge_around_vertex_circulator hec = source->vertex_begin();
	Halfedge_around_vertex_circulator hec_end = hec;

	CGAL_For_all(hec, hec_end)
	{
		if (hec->opposite()->is_border()) break;
	}

	path.clear();

	bool connected = false;

	for (Halfedge_handle he = hec->opposite(); he->vertex() != source; he = he->next())
	{
		path.push_back(he);
		length += he->uvLength();

		if (he->vertex() == target)
		{
			connected = true;
			break;
		}
	}

	if (connected)
	{
		return length;
	}
	else
	{
		return 0.0; //there is no border path from source to target
	}
}



double CGAL_Borders::getPath3DLength(const std::vector<Halfedge_handle>& path)
{
	double pathLength = 0.0;
	for (int i = 0; i < path.size(); i++)
	{
		double edgeLength = path[i]->length(); //edge length from 3D embedding
		pathLength += edgeLength;
	}
	return pathLength;
}



//find the halfedges of the border that the vertex belongs to
//returns the uv length of the polyline border or 0.0 if such a border is not found
double CGAL_Borders::getBorderHalfEdges(const Vertex_handle& vertex, std::vector<Halfedge_handle>& path)
{
	double length = 0.0;

	if (!vertex->is_border())
	{
		return 0.0;
	}

	Halfedge_around_vertex_circulator hec = vertex->vertex_begin();
	Halfedge_around_vertex_circulator hec_end = hec;

	CGAL_For_all(hec, hec_end)
	{
		if (hec->opposite()->is_border()) break;
	}

	path.clear();

	Halfedge_handle he = hec->opposite();

	for (; he->vertex() != vertex; he = he->next())
	{
		path.push_back(he);
		length += he->uvLength();
	}
	path.push_back(he); //add the last one
	length += he->uvLength();

	return length;
}


