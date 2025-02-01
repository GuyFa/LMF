#pragma once


//the functions here assumes that the mesh is manifold.
//meaning it doesn't have singular edges and singular vertices.
//I never tested this code for general polygonal meshes but rather only on triangle meshes.

class CGAL_Borders
{
public:

	CGAL_Borders(Mesh& mesh);
	void printMeshTopologicalStatistics();
	int genus() const;
	int eulerNumber() const;
	int numVertices() const;
	int numFaces() const;
	int numEdges() const;
	int numConnectedComponents() const;
	int numBorders() const;
	int numVerticesInBorder(unsigned int index) const;
	bool isTopologicalDisk() const;
	bool isMeshConnected() const;
	Vertex_handle vertex(unsigned int borderIndex, unsigned int vertexIndex) const;
	double borderCurvature(unsigned int index, bool useTargetMetric) const;
	static double getBorderHalfEdgePath(const Vertex_handle& source, const Vertex_handle& target, std::vector<Halfedge_handle>& path);
	static double getBorderHalfEdges(const Vertex_handle& vertex, std::vector<Halfedge_handle>& path);
	static double getPath3DLength(const std::vector<Halfedge_handle>& path);

private:
	unsigned int computeNumberOfConnectedComponents(Mesh& mesh);
	void traverseConnectedComponent(Facet_handle rootFacet, const double unVisitedTag, const double visitedTag);
	void traverseBorderComponent(Halfedge_handle he, const int unVisitedTag, const int visitedTag);
	void makeBorderHalfedge(Halfedge_iterator& he_border);
	void computeBorders(Mesh& mesh);


private:
	int mV;
	int mF;
	int mE;
	int mC;
	int mB;

	std::vector<std::vector<Vertex_handle> > mBorders;
};
