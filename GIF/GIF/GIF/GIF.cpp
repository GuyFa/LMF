#include "GIF.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include <Eigen/Sparse>
#include "Utils/CompilationAccelerator.h"
#include "Utils/MeshAlgorithms.h"
#include <CGAL/Timer.h>
#include <random>
#include <numeric>


GIF::GIF() {
	std::cout.rdbuf(std::cerr.rdbuf()); //added since cout does not print to output window
}

bool GIF::run(std::string obj_path, std::string uv_mat_path, bool run_IDT, bool allow_MVC_fix, std::string edge_lengths_mat_path, int interior_faces, int boundary_segment_size, double curvature_meta_vertices_rate, double outer_termination_condition_rate, double energy_related_termination_condition)
{
	mNumInternalFaces = interior_faces;
	mSegSize = boundary_segment_size;
	mPercentCurvatureMeta = curvature_meta_vertices_rate;
	mOuterRateTerminationCondition = outer_termination_condition_rate;
	mEnergyRelatedTerminationCondition = energy_related_termination_condition;
	mUseCurvatureMetaVertices = true; // this option adds the vertices with the largest Gaussian curvature to the meta vertices set
	mRunIDT = run_IDT;
	mAllowMVCFix = allow_MVC_fix;
	mDelaunayOutputs.clear();
	mDelaunayOutputs.resize(10, 1);

	//Load mesh
	bool res = loadMesh(obj_path, edge_lengths_mat_path);
	if (!res) {
		cout << "Failed to load mesh" << endl;
		return false;
	}
	//Get borders of mesh
	CGAL_Borders borders(mCgalMesh);
	if (borders.genus() != 0) {
		cout << "The code only supports genus 0 for now" << endl;
		return false;
	}
	if (borders.numBorders() != 1) {
		cout << "The code only supports models with one border" << endl;
		return false;
	}
	//Initialize mesh and arrays, and sent to MATLAB
	res = initialize(borders);
	if (!res) {
		cout << "Failed to initialize mesh" << endl;
		return false;
	}

	//Create KKT, calculate basis, run Newton, and find all UV's
	CGAL::Timer totalTime;
	totalTime.start();
	res = runAlgorithm(borders);
	if (!res) {
		cout << "Failed in algorithm" << endl;
		return false;
	}

	//Update UV's for halfedges in the mesh
	updatHalfedgeUVs();

	//Decide if failed, succeeded or partially succeeded
	getResult();
	if (mPerformFinalGlobalScaling)
	{
		CGAL::Timer globalScalingTime;
		globalScalingTime.start();
		MatlabInterface::GetEngine().EvalToCout("global_scaling");
		globalScalingTime.stop();
		cout << "global scaling time: " << globalScalingTime.time() << endl;
		MatlabGMMDataExchange::GetEngineDenseMatrix("array1", mUVs);
		totalTime.stop();
		updatHalfedgeUVs();
	}
	else
		totalTime.stop();
	
	std::cout << "Total time: " << totalTime.time() << "\n";
	mParameters_matrix(4, 0) = totalTime.time();

	GMMDenseColMatrix uvs(mCgalMesh.size_of_vertices(), 2);
	for_each_vertex(vi, mCgalMesh)
	{
		uvs(vi->index(), 0) = real(vi->uvC());
		uvs(vi->index(), 1) = imag(vi->uvC());
	}

	vector<double> absolute_errors, relative_errors;

	for (auto he = mCgalMesh.edges_begin(); he != mCgalMesh.edges_end(); he++)
		if (!he->is_border_edge()) {
			double first_edge_length = he->targetMetric();
			double second_edge_length = he->opposite()->targetMetric();
			double absolute_error = abs(first_edge_length - second_edge_length);
			double relative_error = absolute_error / ((first_edge_length + second_edge_length) / 2);
			absolute_errors.push_back(absolute_error);
			relative_errors.push_back(relative_error);
		}

	double absolute_erros_sum = 0.0;
	for (double num : absolute_errors) {
		absolute_erros_sum += num;
	}
	double absolute_erros_average = absolute_erros_sum / absolute_errors.size();
	double sum = 0.0;
	for (double num : absolute_errors) {
		sum += (num - absolute_erros_average) * (num - absolute_erros_average);
	}
	double absolute_erros_std = std::sqrt(sum / absolute_errors.size());

	mParameters_matrix(15, 0) = absolute_erros_average;
	mParameters_matrix(16, 0) = absolute_erros_std;


	double relative_erros_sum = 0.0;
	for (double num : relative_errors) {
		relative_erros_sum += num;
	}
	double relative_erros_average = relative_erros_sum / relative_errors.size();
	sum = 0.0;
	for (double num : relative_errors) {
		sum += (num - relative_erros_average) * (num - relative_erros_average);
	}
	double relative_erros_std = std::sqrt(sum / relative_errors.size());

	mParameters_matrix(17, 0) = relative_erros_average;
	mParameters_matrix(18, 0) = relative_erros_std;


	MatlabGMMDataExchange::SetEngineDenseMatrix("uvs_ordered", uvs);
	MatlabGMMDataExchange::SetEngineDenseMatrix("flips_count_before_fix", mFlipsCountBeforeFix);
	MatlabGMMDataExchange::SetEngineDenseMatrix("flips_count_after_fix", mFlipsCountAfterFix);
	MatlabGMMDataExchange::SetEngineDenseMatrix("IDT_outputs", mDelaunayOutputs);
	MatlabGMMDataExchange::SetEngineDenseMatrix("parametersMatrix", mParameters_matrix);
	
	std::string cmd_save_uv = "save('" + uv_mat_path + "', 'uvs_ordered', 'flips_count_before_fix', 'flips_count_after_fix', 'IDT_outputs', 'parametersMatrix');";
	MatlabInterface::GetEngine().Eval(cmd_save_uv.c_str());
	cout << "GIF succeed!" << endl;

	return true;
}


bool GIF::initialize(CGAL_Borders& borders)
{
	UpdateIndexOfHEinSystem();
	getBordersMapAndSetMetaVertices(borders);
	bool hasBorder = mConesAndMetaMap.size() > 0;
	if (!hasBorder) {
		cout << "Model must have border" << endl;
		return false;
	}

	int numOfBorderVertices = 0;
	for (int i = 0; i < borders.numBorders(); i++)
		numOfBorderVertices += borders.numVertices();
	mNumDOF = hasBorder ? mConesAndMetaMap.size() : mConesAndMetaMap.size() - 1;
	mNumScaffoldDOF = min(max((int)(1.1 * (double)mNumDOF), (mNumDOF + 3)), (25 + mNumDOF));
	bool res = TransferBorderFacesAndVerticesToMatlab();
	if (!res) {
		cout << "Failed to transfer to MATLAB" << endl;
		return false;
	}
	return true;
}



bool GIF::loadMesh(std::string obj_path, std::string edge_lengths_mat_path)
{
	// -----load mesh-------
	MeshBuffer m;
	mObjpath = obj_path;
	//MatlabInterface::GetEngine().Eval("GIF.meshName = objLocation;");

	//mSegSize = 1;
	mOperationMode = 1;
	mPerformFinalGlobalScaling = false;

	cout << "Loading mesh information, creating CGAL mesh and transferring mesh to matlab..." << endl;
	bool res = Parser::loadOBJ(mObjpath.c_str(), &mMeshBuffer, &m);
	if (!res) return false;

	// -----create cgal mesh-------
	std::vector<Point_3> pVec;
	pVec.resize(mMeshBuffer.positions.size());
	for (int i = 0; i < (int)pVec.size(); ++i)
		pVec[i] = Point_3(mMeshBuffer.positions[i][0], mMeshBuffer.positions[i][1], mMeshBuffer.positions[i][2]);
	std::vector<unsigned int> faces;
	faces.resize((int)mMeshBuffer.idx_pos.size());
	for (int i = 0; i < (int)faces.size(); ++i)
		faces[i] = mMeshBuffer.idx_pos[i];
	MeshBuilder<HDS, Kernel> meshBuilder(&pVec, (std::vector<int>*)(&faces), NULL);//we assume here that the indices are not negative!
	mCgalMesh.clear();
	mCgalMesh.delegate(meshBuilder);
	mCgalMesh.updateAllGlobalIndices();


	if (mSegSize == -1)
		if (mCgalMesh.size_of_facets() > 500000)
			mSegSize = ceil(double(mCgalMesh.size_of_border_edges()) / 250);
		else
			mSegSize = 10;
		if (mCgalMesh.size_of_border_edges() < 100)
			mSegSize = 1;
	
	if (mNumInternalFaces == -1)
		if (mCgalMesh.size_of_facets() < 8000)
			mNumInternalFaces = 50;
		else if (mCgalMesh.size_of_facets() >= 8000 && mCgalMesh.size_of_facets() < 250000)
			mNumInternalFaces = 300;
		else
			mNumInternalFaces = 500;


	if (edge_lengths_mat_path == "")
	{
		mCgalMesh.updateTargetMetricFrom3dEmbedding();
		cout << "Using original mesh's edge lengths" << endl;
	}
	else
	{
		std::string cmd = "load('" + edge_lengths_mat_path + "')";
		MatlabInterface::GetEngine().EvalToCout(cmd.c_str());
		GMMDenseColMatrix edge_lengths(mCgalMesh.size_of_facets(), 3);
		MatlabGMMDataExchange::GetEngineDenseMatrix("edge_lengths", edge_lengths);

		int j = 0;
		for_each_facet(f, mCgalMesh)
		{
			Halfedge_handle h1 = f->halfedge();
			Halfedge_handle h2 = h1->next();
			Halfedge_handle h3 = h2->next();

			h1->targetMetric() = edge_lengths(j, 2);
			h2->targetMetric() = edge_lengths(j, 0);
			h3->targetMetric() = edge_lengths(j, 1);

			if (h1->is_border_edge())
				h1->opposite()->targetMetric() = edge_lengths(j, 2);
			if (h2->is_border_edge())
				h2->opposite()->targetMetric() = edge_lengths(j, 0);
			if (h3->is_border_edge())
				h3->opposite()->targetMetric() = edge_lengths(j, 1);
			j++;
		}
		cout << "Using edge lengths from NJF" << endl;
	}
	if (mRunIDT)
		for (auto he = mCgalMesh.edges_begin(); he != mCgalMesh.edges_end(); he++)
			if (!he->is_border_edge()) {
				double first_edge_length = he->targetMetric();
				double second_edge_length = he->opposite()->targetMetric();
				double mean_edge_length = (first_edge_length + second_edge_length) / 2.0;
				he->targetMetric() = mean_edge_length;
				he->opposite()->targetMetric() = mean_edge_length;
			}
	if (mRunIDT)
		hyperbolicDelaunayTriangulation(false);
	return true;
}
void GIF::getBordersMapAndSetMetaVertices(CGAL_Borders& borders)
{
	mConesAndMetaMap.clear();
	mIndicesOfMetaVerticesInUVbyGeneralIndex.clear();
	mIndicesOfMetaVerticesInUVbyGeneralIndex.resize(mCgalMesh.size_of_vertices());
	mIndicesOfAllBorderVerticesInUVbyGeneralIndex.clear();
	mIndicesOfAllBorderVerticesInUVbyGeneralIndex.resize(mCgalMesh.size_of_vertices());
	mGeneralIndicesUVbyIndicesOfAllBorderVertices.clear();
	mGeneralIndicesUVbyIndicesOfAllBorderVertices.resize(mCgalMesh.size_of_border_edges());
	mIsMeta.clear();
	mIsMeta.resize(mCgalMesh.size_of_vertices(), false);
	getVerticesThatMustBeMeta(borders);
	mMetaVerticesInBorderByHalfEdge.clear();
	mMetaVerticesInBorderByHalfEdge.resize(borders.numBorders());
	mConesVertices.clear();
	int ind1 = 0;
	//find meta vertices

	int borderSegmentSize = mSegSize;
	std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[0];

	std::vector<Halfedge_handle> borderPath; //halfedges on boundary
	borders.getBorderHalfEdges(borders.vertex(0, 0), borderPath);
	int min_vertex_index = 0;
	for (int j = 0; j < borderPath.size(); j++)
		if (borderPath[j]->vertex()->index() < borderPath[min_vertex_index]->vertex()->index())
			min_vertex_index = j;
	std::rotate(borderPath.begin(), borderPath.begin() + (min_vertex_index + 1), borderPath.end());



	int curvatueSamples = int(mPercentCurvatureMeta * borderPath.size());
	std::vector<double> gaussianCurvatures;
	for (int j = 0; j < borderPath.size(); j++) {

		Halfedge_handle HD = borderPath[j];
		gaussianCurvatures.push_back(abs(HD->vertex()->gaussianCurvature()));
	}
	for (int j = 0; j < curvatueSamples; j++) {
		auto i = std::max_element(gaussianCurvatures.begin(), gaussianCurvatures.end());
		int i1 = std::distance(gaussianCurvatures.begin(), i);
		Halfedge_handle HD = borderPath[i1];
		mIsMeta[HD->vertex()->index()] = true;
		gaussianCurvatures[i1] = 0.0;
	}


	if (borderPath.size() < 3 * mSegSize) {//make sure we have at least 3 meta vertices in this border
		borderSegmentSize = borderPath.size() / 3;
	}

	int verticesInSeg = borderSegmentSize;
	for (int j = borderPath.size() - 1; j >= 0; j--) {

		Halfedge_handle HD = borderPath[j];
		mIndicesOfAllBorderVerticesInUVbyGeneralIndex[HD->vertex()->index()] = ind1;
		mGeneralIndicesUVbyIndicesOfAllBorderVertices[ind1] = HD->vertex()->index();
		mBoundaryVertices.insert(mIndexOfHEinSystem[HD->index()]);
		mBoundaryVerticesMap.push_back(HD->vertex());
		ind1++;
		if (verticesInSeg >= borderSegmentSize || mIsMeta[HD->vertex()->index()]) {//meta vertex
			currBorder.push_back(HD);
			verticesInSeg = 1;
			mIndicesOfMetaVerticesInUVbyGeneralIndex[HD->vertex()->index()] = mConesAndMetaMap.size();
			mConesAndMetaMap.push_back(HD->vertex());
			mConesVertices.insert(mIndexOfHEinSystem[HD->index()]);
			mIsMeta[HD->vertex()->index()] = true;
		}
		else
			verticesInSeg++;
	}
	currBorder.push_back(currBorder[0]);

}

void GIF::hyperbolicDelaunayTriangulation(bool allow_non_planar_edge_flips)
{
	// Retriangulating using Ptolemy flips
	// Assuming edge lengths are consistent on both sides of each edge!
	// Triangle inequality doesn't have to be satisfied
	cout << "Applying hyperbolic Delaunay Triangulation..." << endl;
	const double epsilon = 2e-14;
	CGAL::Timer idtTime;
	idtTime.start();
	int bDomain = 0;
	vector<bool> bMarked(mCgalMesh.size_of_halfedges() / 2, false);
	stack<Halfedge_handle> s;
	mCgalMesh.updateTargetMetricFrom3dEmbedding();
	for (auto he = mCgalMesh.edges_begin(); he != mCgalMesh.edges_end(); he++)
		if (!he->is_border_edge()) {
			s.push(he);
			bMarked[he->index() / 2] = true;
		}
	int nFlips = 0;
	int nUnabledFlips = 0;
	while (!s.empty()) {
		auto he = s.top();
		s.pop();
		bMarked[he->index() / 2] = false;

		double l_mi = he->prev()->targetMetric();
		double l_mj = he->next()->targetMetric();
		double l_kj = he->opposite()->prev()->targetMetric();
		double l_ki = he->opposite()->next()->targetMetric();
		double l_ij = he->targetMetric();

		bool is_triangle_inequality_satisfied = ((l_ki + l_kj - l_ij > 0) && (l_ki - l_kj + l_ij > 0) && (-l_ki + l_kj + l_ij > 0)) && ((l_mi + l_mj - l_ij > 0) && (l_mi - l_mj + l_ij > 0) && (-l_mi + l_mj + l_ij > 0));

		bool local_ideal_Delaunay_condition_holds = ((l_ki * l_ki + l_kj * l_kj - l_ij * l_ij) / (2 * l_ki * l_kj) + (l_mi * l_mi + l_mj * l_mj - l_ij * l_ij) / (2 * l_mi * l_mj)) >= -epsilon;

		if (local_ideal_Delaunay_condition_holds)
			continue;

		// flip
		const double l_km = (l_ki * l_mj + l_kj * l_mi) / l_ij;

		

		he = mCgalMesh.flip_edge(he);

		if (!allow_non_planar_edge_flips)
		{
			//here we check if the mesh remains planar after the edge flip

			if (he->vertex()->index() == he->opposite()->vertex()->index())  // here we check if the new edge connects a vertex to itself
			{
				cout << "irregular mesh operation" << endl;
				he = mCgalMesh.flip_edge(he);
				he->targetMetric() = l_ij;
				he->opposite()->targetMetric() = l_ij;
				nUnabledFlips++;
			}


			// here we check if there exist 2 edges connecting the same 2 vertices, creating a "hat" which is non-planar
			Mesh::Halfedge_around_vertex_circulator h = he->vertex()->vertex_begin();
			const Mesh::Halfedge_around_vertex_circulator hEnd = h;
			int sameEdgesAfterFlip = 0;
			CGAL_For_all(h, hEnd)
			{
				if (he->opposite()->vertex()->index() == h->opposite()->vertex()->index())
					sameEdgesAfterFlip++;
			}
			if (sameEdgesAfterFlip != 1)
			{
				cout << sameEdgesAfterFlip << " edges between 2 vertices!" << endl;
				he = mCgalMesh.flip_edge(he);
				he->targetMetric() = l_ij;
				he->opposite()->targetMetric() = l_ij;
				nUnabledFlips++;

				continue;
			}
		}

		// after the flip the edge is in the direction of the triangle
		he->targetMetric() = l_km;
		he->opposite()->targetMetric() = l_km;

		bool local_ideal_Delaunay_condition_holds_after_flip = ((l_ki * l_ki + l_mi * l_mi - l_km * l_km) / (2 * l_ki * l_mi) + (l_kj * l_kj + l_mj * l_mj - l_km * l_km) / (2 * l_kj * l_mj)) >= 0;
		assert(local_ideal_Delaunay_condition_holds_after_flip);
		bool is_triangle_inequality_satisfied_1 = ((l_mj + l_kj - l_km > 0) && (l_mj - l_kj + l_km > 0) && (-l_mj + l_kj + l_km > 0)) && ((l_mi + l_ki - l_km > 0) && (l_mi - l_ki + l_km > 0) && (-l_mi + l_ki + l_km > 0));
		assert(is_triangle_inequality_satisfied_1);
		nFlips++;

		// add four neighbors to stack
		vector<Halfedge_handle> ne;
		for (int i = 0; i < 2; ++i) {

			ne.push_back(he->next());
			ne.push_back(he->prev());
			he = he->opposite();
		}
		for (auto he : ne)
		{
			if (!he->is_border_edge() && !bMarked[he->index() / 2])
			{
				s.push(he);
				bMarked[he->index() / 2] = true;
			}
		}
	}
	idtTime.stop();
	cout << "The number of flips IDT performed is " << nFlips << endl;
	cout << "Intrinsic Delaunay Triangulation took " << idtTime.time() << endl;
	int negativeCotEdgesCount = 0, internalEdgesCount = 0, nonIdealCriteriaSatisfyingEdgesCount = 0;
	double minCotValue = 1.0, negativeCotEdgesPercent = 0.0;
	bool areAllCotPositive = true, areIdealCriteriasHold = true;

	mCgalMesh.updateAllGlobalIndices();
	
	for (auto heIt1 = mCgalMesh.edges_begin(); heIt1 != mCgalMesh.edges_end(); heIt1++)
	{
		if (heIt1->is_border_edge()) // checks both current halfedge or its opposite, for border edge
			continue;
		internalEdgesCount++;
		if ((((heIt1->cot(true) + heIt1->opposite()->cot(true))) / 2) < minCotValue)
			minCotValue = ((heIt1->cot(true) + heIt1->opposite()->cot(true))) / 2;
		if (((heIt1->cot(true) + heIt1->opposite()->cot(true)) < -epsilon))
		{
			negativeCotEdgesCount++;
			areAllCotPositive = false;
		}
		double l_mi = heIt1->prev()->targetMetric();
		double l_mj = heIt1->next()->targetMetric();
		double l_kj = heIt1->opposite()->prev()->targetMetric();
		double l_ki = heIt1->opposite()->next()->targetMetric();
		double l_ij = heIt1->targetMetric();
		bool local_ideal_Delaunay_condition_holds = ((l_ki * l_ki + l_kj * l_kj - l_ij * l_ij) / (2 * l_ki * l_kj) + (l_mi * l_mi + l_mj * l_mj - l_ij * l_ij) / (2 * l_mi * l_mj)) >= -epsilon;
		if (!local_ideal_Delaunay_condition_holds)
		{
			nonIdealCriteriaSatisfyingEdgesCount++;
			areIdealCriteriasHold = false;
		}
	}
	negativeCotEdgesPercent = double(negativeCotEdgesCount) / double(internalEdgesCount);
	cout << "are all cot positive " << areAllCotPositive << endl;
	cout << "do all ideal Delaunay criterias hold " << areIdealCriteriasHold << endl;
	cout << "minimum cot weight is " << minCotValue << endl;
	cout << "negative cot edges count is " << negativeCotEdgesCount << endl;
	cout << "non ideal criteria satisfying edges count is " << nonIdealCriteriaSatisfyingEdgesCount << endl;
	cout << "internal edges count is " << internalEdgesCount << endl;
	cout << "negative cot edges rate is " << negativeCotEdgesPercent << endl;
	bool isSimpleGraph = true;
	for_each_vertex(v1, mCgalMesh)
	{
		Mesh::Halfedge_around_vertex_circulator h = v1->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;
		std::vector<int> neighbors;
		CGAL_For_all(h, hEnd)
		{
			neighbors.push_back(h->opposite()->vertex()->index());
			if (h->vertex()->index() == h->opposite()->vertex()->index())
				isSimpleGraph = false;
		}
		const int originalNeighborsCount = neighbors.size();
		sort(neighbors.begin(), neighbors.end());
		neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
		if (originalNeighborsCount != neighbors.size())
			isSimpleGraph = false;
	}

	int triangle_inequality_not_satisfied_count = 0;
	auto facetIter = mCgalMesh.facets_begin();
	while (facetIter != mCgalMesh.facets_end())
	{
		Halfedge_handle he[3];
		facetIter->getHalfedges(he);
		double l_mj = he[0]->targetMetric();
		double l_mi = he[1]->targetMetric();
		double l_ij = he[2]->targetMetric();
		bool is_triangle_inequality_satisfied = ((l_mj + l_mi > l_ij) && (l_mj + l_ij > l_mi) && (l_mi + l_ij > l_mj));
		if (!is_triangle_inequality_satisfied)
		{
			triangle_inequality_not_satisfied_count++;
			cout << "l_mj = " << l_mj << endl;
			cout << "l_mi = " << l_mi << endl;
			cout << "l_ij = " << l_ij << endl;
			cout << "problem! triangle inequality is not satisfied" << endl;
		}
		facetIter++;
	}
	cout << "is simple graph " << isSimpleGraph << endl;
	mDelaunayOutputs(0, 0) = nFlips;
	mDelaunayOutputs(1, 0) = nUnabledFlips;
	mDelaunayOutputs(2, 0) = areAllCotPositive;
	mDelaunayOutputs(3, 0) = minCotValue;
	mDelaunayOutputs(4, 0) = negativeCotEdgesCount;
	mDelaunayOutputs(5, 0) = internalEdgesCount;
	mDelaunayOutputs(6, 0) = negativeCotEdgesPercent;
	mDelaunayOutputs(7, 0) = triangle_inequality_not_satisfied_count;
	mDelaunayOutputs(8, 0) = idtTime.time();
	mDelaunayOutputs(9, 0) = isSimpleGraph;
}


void GIF::getScaffoldBordersMapAndSetMetaVertices(CGAL_Borders& borders, int numVertices)
{
	mIsMetaScaffold.clear();
	mIsMetaScaffold.resize(numVertices, false);
	mScaffoldConesAndMetaMap.clear();
	mScaffoldMetaVerticesInBorderByHalfEdge.clear();
	mScaffoldMetaVerticesInBorderByHalfEdge.resize(borders.numBorders());

	double first_border_vertex_x = borders.vertex(0, 0)->point().x();
	double first_border_vertex_y = borders.vertex(0, 0)->point().y();
	double second_border_vertex_x = borders.vertex(1, 0)->point().x();
	double second_border_vertex_y = borders.vertex(1, 0)->point().y();
	double first_border_radius_squared = first_border_vertex_x * first_border_vertex_x + first_border_vertex_y * first_border_vertex_y;
	double second_border_radius_squared = second_border_vertex_x * second_border_vertex_x + second_border_vertex_y * second_border_vertex_y;
	int index_of_inner_border = (first_border_radius_squared < second_border_radius_squared) ? 0 : 1;
	int index_of_outer_border = 1 - index_of_inner_border;


	std::vector<Halfedge_handle>& curr_outer_border = mScaffoldMetaVerticesInBorderByHalfEdge[index_of_outer_border];
	std::vector<Halfedge_handle> outer_border_path; //halfedges on boundary
	borders.getBorderHalfEdges(borders.vertex(index_of_outer_border, 0), outer_border_path);
	std::rotate(outer_border_path.begin(), outer_border_path.begin() + (outer_border_path[0]->vertex()->index() + 1), outer_border_path.end());
	for (int j = outer_border_path.size() - 1; j >= 0; j--)
	{
		Halfedge_handle HD = outer_border_path[j];
		curr_outer_border.push_back(HD);
		mScaffoldConesAndMetaMap.push_back(HD->vertex());//update borders map
		mIsMetaScaffold[HD->vertex()->index()] = true;
	}
	curr_outer_border.push_back(curr_outer_border[0]);

	std::vector<Halfedge_handle>& curr_inner_border = mScaffoldMetaVerticesInBorderByHalfEdge[index_of_inner_border];
	std::vector<Halfedge_handle> inner_border_path; //halfedges on boundary
	borders.getBorderHalfEdges(borders.vertex(index_of_inner_border, 0), inner_border_path);
	int rotateinnerBorder = inner_border_path[inner_border_path.size() - 1]->vertex()->index() - borders.numVerticesInBorder(index_of_outer_border);
	std::rotate(inner_border_path.rbegin(), inner_border_path.rbegin() + rotateinnerBorder, inner_border_path.rend());

	for (int j = inner_border_path.size() - 1; j >= 0; j--)
	{
		Halfedge_handle HD = inner_border_path[j];
		curr_inner_border.push_back(HD);
		mScaffoldConesAndMetaMap.push_back(HD->vertex());//update borders map
		mIsMetaScaffold[HD->vertex()->index()] = true;
	}
	curr_inner_border.push_back(curr_inner_border[0]);
}


void GIF::getVerticesThatMustBeMeta(CGAL_Borders& borders)
{
	// mark vertices on bridges
	GMMDenseColMatrix isVertexOnBridge(mCgalMesh.size_of_vertices(), 1); //indicate if a vertex in the mesh is on a bridge
	gmm::clear(isVertexOnBridge);

	int numBorders = borders.numBorders();
	std::vector<Vertex_handle> firstBridgeVertexInBorder(numBorders, NULL);


	for (int i = 0; i < borders.numBorders(); i++) {
		for (int j = 0; j < borders.numVerticesInBorder(i); j++) {
			Vertex_handle v = borders.vertex(i, j);

			//if a vertex is on cut and on border, is should be meta vertex
			if (v->onCut()) {
				mIsMeta[v->index()] = true;
				continue;
			}

			//find bridges (edges whose deletion disconnects the graph) - and set the 2 vertices on the bridge to be meta vertices
			Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
			const Mesh::Halfedge_around_vertex_circulator hEnd = h;

			double sumWeights = 0.0;
			CGAL_For_all(h, hEnd)
			{
				if (!h->is_border_edge()) {
					Vertex_handle v1 = h->vertex();
					Vertex_handle v2 = h->opposite()->vertex();
					if (v1->is_border() && v2->is_border())
					{
						isVertexOnBridge(v1->index(), 0) = 1;
						isVertexOnBridge(v2->index(), 0) = 1;

						//add bridge vertices to meta
						mIsMeta[v1->index()] = true;
						mIsMeta[v2->index()] = true;

						if (firstBridgeVertexInBorder[i] == NULL)
							firstBridgeVertexInBorder[i] = v;
					}
				}
			}
		}
	}

	//select a vertex between every 2 bridge vertices.
	for (int i = 0; i < borders.numBorders(); i++) {

		if (firstBridgeVertexInBorder[i] == NULL)
			continue;

		Halfedge_around_vertex_circulator hec = firstBridgeVertexInBorder[i]->vertex_begin();
		Halfedge_around_vertex_circulator hec_end = hec;

		CGAL_For_all(hec, hec_end)
		{
			if (hec->is_border()) break;
		}
		Halfedge_handle he = hec;

		do
		{
			Halfedge_handle he_next = he->next();
			while (!mIsMeta[he_next->vertex()->index()])
				he_next = he_next->next();
			Halfedge_handle he2 = he_next;
			while (he->vertex() != he2->opposite()->vertex() && he->next()->vertex() != he2->opposite()->vertex()) {
				he = he->next();
				he2 = he2->prev();
			}
			mIsMeta[he2->opposite()->vertex()->index()] = true;

			he = he_next;
		} while ((he->vertex() != firstBridgeVertexInBorder[i]));
	}
}

void GIF::UpdateIndexOfHEinSystem()
{
	mIndexOfHEinSystem.clear();
	mIndexOfHEinSystem.resize(mCgalMesh.size_of_halfedges());
	mIndexOfVertexByIndexOfHEinSystem.clear();
	mIndexOfVertexByIndexOfHEinSystem.resize(mCgalMesh.size_of_vertices());
	auto heIt = mCgalMesh.halfedges_begin();
	int index1, index2;
	std::vector<bool> visit;
	visit.resize(mCgalMesh.size_of_halfedges());
	for (int i = 0; i < (int)visit.size(); ++i)
		visit[i] = false;
	int num = 0;
	bool flag = false;
	while (heIt != mCgalMesh.halfedges_end())
	{
		if (visit[heIt->index()])
		{
			heIt++;
			continue;
		}

		if (!heIt->vertex()->onCut())//if the vertex is not on the cut, then all the halfedges are mapped to the same variable
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			auto start = oneRing;
			index1 = oneRing->index();
			mIndexOfHEinSystem[index1] = num;
			mIndexOfVertexByIndexOfHEinSystem[num] = heIt->vertex()->index();
			visit[index1] = true;
			oneRing++;
			while (oneRing != start)
			{
				index2 = oneRing->index();
				mIndexOfHEinSystem[index2] = num;
				visit[index2] = true;
				oneRing++;
			}
			num++;
		}
		else//if it is on the cut, split the vertex
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			while ((!oneRing->isCut()) && (!oneRing->opposite()->is_border()))
				oneRing++;
			auto start = oneRing;
			index1 = start->index();
			mIndexOfHEinSystem[index1] = num;
			visit[index1] = true;
			oneRing++;

			while (1)	//its o.k, the break condition have to be true
			{
				bool endFlag = false;
				while (!oneRing->isCut())
				{
					flag = true;
					index2 = oneRing->index();
					mIndexOfHEinSystem[index2] = num;
					visit[index2] = true;
					if (oneRing->is_border())
					{
						if (oneRing == start)
							endFlag = true;
						oneRing++;
						break;
					}
					oneRing++;
				}
				num++;
				if ((oneRing == start) || (endFlag))
					break;
				index1 = oneRing->index();
				mIndexOfHEinSystem[index1] = num;
				visit[index1] = true;
				oneRing++;
			}
		}
		heIt++;
	}
	mSizeOfSystemVar = num;
	GMMDenseColMatrix numOfsystemVars(1, 1);
	numOfsystemVars(0, 0) = mSizeOfSystemVar;
	MatlabGMMDataExchange::SetEngineDenseMatrix("numOfsystemVars", numOfsystemVars);
	MatlabInterface::GetEngine().EvalToCout("GIF.sysVars = numOfsystemVars(1, 1);");
}

void GIF::UpdateIndexOfHEinScaffoldSystem(Mesh& mesh1)
{

	mIndexOfHEinScaffoldSystem.clear();
	mIndexOfHEinScaffoldSystem.resize(mesh1.size_of_halfedges());
	auto heIt = mesh1.halfedges_begin();
	int index1, index2;
	std::vector<bool> visit;
	visit.resize(mesh1.size_of_halfedges());
	for (int i = 0; i < (int)visit.size(); ++i)
		visit[i] = false;
	int num = 0;
	bool flag = false;
	while (heIt != mesh1.halfedges_end())
	{
		if (visit[heIt->index()])
		{
			heIt++;
			continue;
		}

		if (!heIt->vertex()->onCut())//if the vertex is not on the cut, then all the halfedges are mapped to the same variable
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			auto start = oneRing;
			index1 = oneRing->index();
			mIndexOfHEinScaffoldSystem[index1] = num;
			visit[index1] = true;
			oneRing++;
			while (oneRing != start)
			{
				index2 = oneRing->index();
				mIndexOfHEinScaffoldSystem[index2] = num;
				visit[index2] = true;
				oneRing++;
			}
			num++;
		}
		else//if it is on the cut, split the vertex
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			while ((!oneRing->isCut()) && (!oneRing->opposite()->is_border()))
				oneRing++;
			auto start = oneRing;
			index1 = start->index();
			mIndexOfHEinScaffoldSystem[index1] = num;
			visit[index1] = true;
			oneRing++;

			while (1)	//its o.k, the break condition have to be true
			{
				bool endFlag = false;
				while (!oneRing->isCut())
				{
					flag = true;
					index2 = oneRing->index();
					mIndexOfHEinScaffoldSystem[index2] = num;
					visit[index2] = true;
					if (oneRing->is_border())
					{
						if (oneRing == start)
							endFlag = true;
						oneRing++;
						break;
					}
					oneRing++;
				}
				num++;
				if ((oneRing == start) || (endFlag))
					break;
				index1 = oneRing->index();
				mIndexOfHEinScaffoldSystem[index1] = num;
				visit[index1] = true;
				oneRing++;
			}
		}
		heIt++;
	}
	mSizeOfScaffoldSystemVar = num;
}

bool GIF::TransferBorderFacesAndVerticesToMatlab()
{
	//******* Tranfer duplicated Vertices by halfedges************
	mParameters_matrix.clear();
	mParameters_matrix.resize(20, 1);
	mParameters_matrix(0, 0) = mCgalMesh.size_of_vertices();
	mParameters_matrix(1, 0) = mCgalMesh.size_of_facets();
	mParameters_matrix(2, 0) = mCgalMesh.size_of_border_edges();
	mParameters_matrix(3, 0) = mNumDOF;
	mParameters_matrix(5, 0) = mSegSize;
	GMMDenseColMatrix verticesByHalfedges(mSizeOfSystemVar, 3);
	auto heIt = mCgalMesh.halfedges_begin();
	while (heIt != mCgalMesh.halfedges_end())
	{
		int index = mIndexOfHEinSystem[heIt->index()];
		auto p = heIt->vertex()->point();
		verticesByHalfedges(index, 0) = p[0];
		verticesByHalfedges(index, 1) = p[1];
		verticesByHalfedges(index, 2) = p[2];
		heIt++;
	}

	GMMDenseColMatrix facesIndices(mCgalMesh.size_of_facets(), 3);
	std::vector<int> boundaryFacesIndicesForCheck;
	std::vector<int> internalFacesIndicesForCheck;
	auto facetIter = mCgalMesh.facets_begin();
	int i = 0, boundaryFacesIndex = 0, internalFacesIndex = 0;
	while (facetIter != mCgalMesh.facets_end())
	{
		Halfedge_handle he[3];
		facetIter->getHalfedges(he);

		Vertex_handle vertices1[3];
		facetIter->getVertices(vertices1);

		facesIndices(i, 0) = mIndexOfHEinSystem[he[0]->index()] + 1;
		facesIndices(i, 1) = mIndexOfHEinSystem[he[1]->index()] + 1;
		facesIndices(i, 2) = mIndexOfHEinSystem[he[2]->index()] + 1;
		if (vertices1[0]->is_border() || vertices1[1]->is_border() || vertices1[2]->is_border())
			boundaryFacesIndicesForCheck.push_back(i + 1); // boundary faces
		else
			internalFacesIndicesForCheck.push_back(i + 1); // internal faces

		facetIter++;
		i++;
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.F1", facesIndices);

	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.V", verticesByHalfedges);

	//********* Transfer border faces corresponding to Halfedges and vector field****************
	std::unordered_set<Facet_handle> facesNearCones;//use hash to prevent duplications
	std::unordered_set<Facet_handle> facesNearCones1;

	for (int i = 0; i < mConesAndMetaMap.size(); i++) {
		Vertex_handle v = mConesAndMetaMap[i];
		int ind = v->index();
		Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;
		CGAL_For_all(h, hEnd)
		{
			if (!h->is_border()) {
				facesNearCones.insert(h->face());
				facesNearCones1.insert(h->face());
			}
		}
	}

	int numOuterFacesAlready = facesNearCones.size();
	if (mNumInternalFaces > 0)
	{
		GMMDenseColMatrix numFaces(1, 1);
		numFaces(0, 0) = mCgalMesh.size_of_facets();
		MatlabGMMDataExchange::SetEngineDenseMatrix("numFaces", numFaces);
		MatlabInterface::GetEngine().EvalToCout("rng(1);");
		MatlabInterface::GetEngine().EvalToCout("randomTrianglesIndices = (randperm(numFaces(1, 1)) - 1).';");
		GMMDenseColMatrix randomTriangles;
		MatlabGMMDataExchange::GetEngineDenseMatrix("randomTrianglesIndices", randomTriangles);
		int i = 0;
		while (true)
		{
			facesNearCones.insert(mCgalMesh.face(randomTriangles(i, 0)));
			i++;
			if ((facesNearCones.size() - numOuterFacesAlready) == mNumInternalFaces || facesNearCones.size() == mCgalMesh.size_of_facets())
				break;

		}
	}

	mNearBoundaryVertices.clear();
	std::unordered_set<Facet_handle> facesNearCones2;
	for (int i = 0; i < mBoundaryVerticesMap.size(); i++) {
		Vertex_handle v = mBoundaryVerticesMap[i];
		int ind = v->index();
		Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;
		CGAL_For_all(h, hEnd)
		{
			if (!h->is_border()) {
				facesNearCones2.insert(h->face());
			}
		}
	}
	unordered_set<Facet_handle>::const_iterator itr;
	for (itr = facesNearCones2.begin(); itr != facesNearCones2.end(); ++itr) {
		Facet_handle f = (*itr);
		Halfedge_handle he[3];
		f->getHalfedges(he);

		if (!mBoundaryVertices.count(mIndexOfHEinSystem[he[0]->index()]))
			mNearBoundaryVertices.insert(mIndexOfHEinSystem[he[0]->index()]);
		if (!mBoundaryVertices.count(mIndexOfHEinSystem[he[1]->index()]))
			mNearBoundaryVertices.insert(mIndexOfHEinSystem[he[1]->index()]);
		if (!mBoundaryVertices.count(mIndexOfHEinSystem[he[2]->index()]))
			mNearBoundaryVertices.insert(mIndexOfHEinSystem[he[2]->index()]);

	}
	unordered_set<int>::const_iterator itr1;
	for (itr1 = mNearBoundaryVertices.begin(); itr1 != mNearBoundaryVertices.end(); itr1++)
	{
		mNearBoundaryVerticesVec.push_back((*itr1));
	}

	int numNearConesFaces = facesNearCones.size();
	GMMDenseColMatrix facesByHalfedges(numNearConesFaces, 3);
	GMMDenseColMatrix internalFacesIndicator(numNearConesFaces, 1);

	GMMDenseColMatrix nearConesFaces(numNearConesFaces, 1);

	mConesAndNearConesMapOfRowsInKKT.clear();
	mNearConesVertices.clear();
	i = 0;
	for (itr = facesNearCones.begin(); itr != facesNearCones.end(); ++itr) {
		Facet_handle f = (*itr);
		Halfedge_handle he[3];
		f->getHalfedges(he);
		facesByHalfedges(i, 0) = mIndexOfHEinSystem[he[0]->index()] + 1;//+1 for matlab
		facesByHalfedges(i, 1) = mIndexOfHEinSystem[he[1]->index()] + 1;
		facesByHalfedges(i, 2) = mIndexOfHEinSystem[he[2]->index()] + 1;
		internalFacesIndicator(i, 0) = facesNearCones1.count(f) ? 0 : 1;

		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[0]->index()]);
		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[1]->index()]);
		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[2]->index()]);

		if (!mConesVertices.count(mIndexOfHEinSystem[he[0]->index()]))
			mNearConesVertices.insert(mIndexOfHEinSystem[he[0]->index()]);
		if (!mConesVertices.count(mIndexOfHEinSystem[he[1]->index()]))
			mNearConesVertices.insert(mIndexOfHEinSystem[he[1]->index()]);
		if (!mConesVertices.count(mIndexOfHEinSystem[he[2]->index()]))
			mNearConesVertices.insert(mIndexOfHEinSystem[he[2]->index()]);
		++i;
	}
	mNearConesVerticesVec.clear();
	for (itr1 = mNearConesVertices.begin(); itr1 != mNearConesVertices.end(); itr1++)
	{
		mNearConesVerticesVec.push_back((*itr1));
	}
	mNonBoundaryVerticesNearConesVerticesVec.clear();
	mNonMetaBoundaryVerticesNearConesVerticesVec.clear();
	for (itr1 = mNearConesVertices.begin(); itr1 != mNearConesVertices.end(); itr1++)
	{
		if (!mBoundaryVertices.count((*itr1)))
			mNonBoundaryVerticesNearConesVerticesVec.push_back((*itr1));
		else
			mNonMetaBoundaryVerticesNearConesVerticesVec.push_back((*itr1));
	}
	mMetaIndicesVec.clear();
	mIsNearMeta.clear();
	mIsNearMeta.resize(mCgalMesh.size_of_vertices(), false);
	mMetaIndicesVec.resize(mCgalMesh.size_of_vertices(), -1);
	int num = 0;
	for (int i = 0; i < mConesAndMetaMap.size(); i++)
	{
		mMetaIndicesVec[mConesAndMetaMap[i]->index()] = num;
		num++;
	}
	for (int i = 0; i < mNearConesVerticesVec.size(); i++)
	{
		mIsNearMeta[mNearConesVerticesVec[i]] = true;
	}
	mInternalVerticesIndices.clear();
	mInternalVerticesIndices.resize(mCgalMesh.size_of_vertices());
	num = 0;
	for_each_vertex(vi, mCgalMesh)
	{
		if (!vi->is_border())
		{
			mInternalVerticesIndices[vi->index()] = num;
			num++;
		}
	}
	num = 0;
	mBoundaryVerticesIndices.clear();
	mBoundaryVerticesIndices.resize(mCgalMesh.size_of_vertices());
	for (int i = 0; i < mBoundaryVerticesMap.size(); i++)
	{
		mBoundaryVerticesIndices[mBoundaryVerticesMap[i]->index()] = num;
		num++;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.F", facesByHalfedges);
	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.internalFacesIndicator", internalFacesIndicator);
	GMMDenseColMatrix borderEdges(mBoundaryVerticesMap.size(), 3);
	for (int i = 0; i < mBoundaryVerticesMap.size(); i++) {
		Vertex_handle v = mBoundaryVerticesMap[i];
		borderEdges(i, 0) = v->point().x();
		borderEdges(i, 1) = v->point().y();
		borderEdges(i, 2) = v->point().z();
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.borderEdges", borderEdges);
	mParameters_matrix(6, 0) = facesNearCones1.size();
	mParameters_matrix(7, 0) = facesNearCones.size() - facesNearCones1.size();

	return true;
}

bool GIF::TransferScaffoldBorderFacesAndVerticesToMatlab(Mesh& mesh1)
{

	//******* Tranfer duplicated Vertices by halfedges************
	GMMDenseColMatrix verticesByHalfedges(mesh1.size_of_vertices(), 3);

	auto heIt = mesh1.halfedges_begin();
	while (heIt != mesh1.halfedges_end())
	{
		int index = mIndexOfHEinScaffoldSystem[heIt->index()];
		auto p = heIt->vertex()->point();
		verticesByHalfedges(index, 0) = p[0];
		verticesByHalfedges(index, 1) = p[1];
		verticesByHalfedges(index, 2) = p[2];
		heIt++;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.SV ", verticesByHalfedges);

	//********* Transfer border faces corresponding to Halfedges and vector field****************
	std::unordered_set<Facet_handle> facesNearCones;//use hash to prevent duplications

	for (int i = 0; i < mScaffoldConesAndMetaMap.size(); i++) {
		Vertex_handle v = mScaffoldConesAndMetaMap[i];
		Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;
		CGAL_For_all(h, hEnd)
		{
			if (!h->is_border()) {
				facesNearCones.insert(h->face());
			}
		}
	}
	int numNearConesFaces = facesNearCones.size();
	GMMDenseColMatrix facesByHalfedges(numNearConesFaces, 3);


	unordered_set<Facet_handle>::const_iterator itr;
	int i = 0;
	for (itr = facesNearCones.begin(); itr != facesNearCones.end(); ++itr) {
		Facet_handle f = (*itr);
		Halfedge_handle he[3];
		f->getHalfedges(he);
		facesByHalfedges(i, 0) = mIndexOfHEinScaffoldSystem[he[0]->index()] + 1;//+1 for matlab
		facesByHalfedges(i, 1) = mIndexOfHEinScaffoldSystem[he[1]->index()] + 1;
		facesByHalfedges(i, 2) = mIndexOfHEinScaffoldSystem[he[2]->index()] + 1;
		++i;
	}


	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.SF", facesByHalfedges);

	return true;
}


bool GIF::runAlgorithm(CGAL_Borders& borders)
{
	cout << "Starting algorithm" << endl;
	//*** run algorithm to calculate UV's ***
	GMMDenseColMatrix RHS;
	int conesConstraintsStartRow;

	cout << "mSegSize = " << mSegSize << endl;
	cout << "mNumInternalFaces = " << mNumInternalFaces << endl;

	CGAL::Timer t1;
	t1.start();

	bool res = true;
	constructKKTmatrix(borders, conesConstraintsStartRow);
	constructHarmonicBasisAndSendToMATLAB(conesConstraintsStartRow);
	bool tutteSuccess = getTutteInitialValue(borders);
	GMMDenseColMatrix uvsOfBoundary(mNumDOF, 2);
	double outerRadiusOfScaffold1 = 0.0;
	MatlabInterface::GetEngine().EvalToCout("computeOptimizationRescaling");
	MatlabGMMDataExchange::GetEngineDenseMatrix("x10", uvsOfBoundary);

	mInitialUVsOfCones.clear();
	mInitialUVsOfCones.resize(mNumDOF);
	for (int i = 0; i < mNumDOF; i++)
	{
		mInitialUVsOfCones[i] = Complex(uvsOfBoundary(i, 0), uvsOfBoundary(i, 1));

	}
	Mesh initialScaffold;
	MeshAlgorithms::buildScaffoldWithFreeOuterBoundary(mInitialUVsOfCones, false, initialScaffold, (mNumScaffoldDOF - mNumDOF), 10);

	CGAL_Borders borders1(initialScaffold);



	Eigen::SparseMatrix<double, Eigen::RowMajor> scaffoldKKtEigen1, scaffoldKKtEigen2;

	getScaffoldBordersMapAndSetMetaVertices(borders1, initialScaffold.size_of_vertices());
	UpdateIndexOfHEinScaffoldSystem(initialScaffold);
	TransferScaffoldBorderFacesAndVerticesToMatlab(initialScaffold);


	res = constructScaffoldHarmonicBasisAndSendToMATLAB();
	if (!res) {
		cout << "Failed to construct Harmonic Basis for the scaffold" << endl;
		return false;
	}
	int numVerticesScaffold = initialScaffold.size_of_vertices();
	GMMDenseColMatrix scaffoldBoundaryVertices(2 * numVerticesScaffold, 1);
	for (int i = 0; i < numVerticesScaffold; i++)
	{
		scaffoldBoundaryVertices(i, 0) = mScaffoldConesAndMetaMap[i]->point().x();
		scaffoldBoundaryVertices(i + numVerticesScaffold, 0) = mScaffoldConesAndMetaMap[i]->point().y();
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldBoundaryVertices", scaffoldBoundaryVertices);
	MatlabInterface::GetEngine().EvalToCout("[ GIF.ScaffoldUVonCones, GIF.FixedIndices, GIF.FixedValues ] = fixFirstCone( scaffoldBoundaryVertices );");


	bool NewtonSucsess;
	NewtonSucsess = runNewtonWithRemeshedScaffold(conesConstraintsStartRow, RHS);

	res = getAllUVsByRHS(RHS);
	if (!res) {
		cout << "Failed to Calculate UVs" << endl;
		return false;
	}
	t1.stop();
	cout << "Net algorithm without further fixes took " << t1.time() << endl;

	return true;
}



void GIF::constructKKTmatrix(CGAL_Borders& borders, int& conesConstraintsStartRow)
{
	std::vector<Eigen::Triplet<double>> tripletListValues;

	tripletListValues.reserve(6 * mCgalMesh.size_of_vertices() + mNumDOF);

	FillLaplacianInKKT(tripletListValues);

	conesConstraintsStartRow = mSizeOfSystemVar;
	mHarmonicBasisInAllV.resize(mSizeOfSystemVar, mNumDOF);
	FillMetaVerticesConstraintsInKKT(borders, tripletListValues, conesConstraintsStartRow);

	mSizeOfMatrix = conesConstraintsStartRow + mNumDOF;

	FillConesConstraintsInKKT(tripletListValues, conesConstraintsStartRow);
	SetElementsForPARDISO(tripletListValues, conesConstraintsStartRow);

	mKKtEigen.resize(mSizeOfMatrix, mSizeOfMatrix);
	mKKtEigen.setFromTriplets(tripletListValues.begin(), tripletListValues.end());
	mKKtEigen.makeCompressed();
}


void GIF::FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	//fill the matrix by going on the triangles one by one
	for_each_facet(f, mCgalMesh)
	{
		Halfedge_handle h1 = f->halfedge();
		Halfedge_handle h2 = h1->next();
		Halfedge_handle h3 = h2->next();
		
		double term1 = -0.5* h1->cot(true);
		double term2 = -0.5* h2->cot(true);
		double term3 = -0.5* h3->cot(true);

		//v1
		updateTermInLaplacianByHalfEdge(h3, term3, tripletListValues); //v1->v2 and v2->v1

		//v2
		updateTermInLaplacianByHalfEdge(h1, term1, tripletListValues); //v2->v3 and v3->v2

		//v3
		updateTermInLaplacianByHalfEdge(h2, term2, tripletListValues); //v3->v1 and v1->v3
	}
}


void GIF::updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	const int i = mIndexOfHEinSystem[h->prev()->index()];
	const int j = mIndexOfHEinSystem[h->index()];

	if (term != 0) {
		updateTermInLaplacianByIndices(i, j, term, tripletListValues);
	}
}


void GIF::updateTermInLaplacianByIndices(int i, int j, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	tripletListValues.push_back(Eigen::Triplet<double>(i, i, -term));
	tripletListValues.push_back(Eigen::Triplet<double>(j, j, -term));
	// only store upper triangular part
	if (j >= i) { //v1->v2
		tripletListValues.push_back(Eigen::Triplet<double>(i, j, term));
	}
	else {//v2->v1
		tripletListValues.push_back(Eigen::Triplet<double>(j, i, term));
	}
}


void GIF::FillMetaVerticesConstraintsInKKT(CGAL_Borders& borders, std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT)
{
	//calculate weights to meta vertices
	std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[0];
	for (int i = 0; i < currBorder.size() - 1; i++) {
		Halfedge_handle metaVertexHD = currBorder[i];
		Halfedge_handle nextMetaVertexHD = currBorder[i + 1];

		//find meta Edge Length
		Halfedge_handle h = metaVertexHD;
		double metaEdgeLength = 0.0;
		while (h != nextMetaVertexHD) {
			metaEdgeLength += h->length();
			h = h->prev();
		}

		//update meta vertices and t
		h = metaVertexHD->prev();
		double lengthFromFirstMeta = 0.0;
		while (h != nextMetaVertexHD) {
			const int firstMetaVertex = mIndexOfHEinSystem[metaVertexHD->index()];
			const int secondMetaVertex = mIndexOfHEinSystem[nextMetaVertexHD->next()->opposite()->index()];//if the vertex is on cut then this gives the correct halfedge, and if not that it is the same.
			const int curVertex = mIndexOfHEinSystem[h->index()];
			lengthFromFirstMeta += h->next()->length();
			const double t = (1 - lengthFromFirstMeta / metaEdgeLength);

			// if we use Eigen, then we only fill upper triangular part, then only A^t is relevant
			tripletListValues.push_back(Eigen::Triplet<double>(firstMetaVertex, rowInKKT, t));
			tripletListValues.push_back(Eigen::Triplet<double>(secondMetaVertex, rowInKKT, 1 - t));
			tripletListValues.push_back(Eigen::Triplet<double>(curVertex, rowInKKT, -1));
			if (mIsNearMeta[curVertex])
			{
				mHarmonicBasisInAllV(curVertex, mMetaIndicesVec[mIndexOfVertexByIndexOfHEinSystem[firstMetaVertex]]) = t;
				mHarmonicBasisInAllV(curVertex, mMetaIndicesVec[mIndexOfVertexByIndexOfHEinSystem[secondMetaVertex]]) = 1 - t;
			}


			++rowInKKT;
			h = h->prev();
		}
	}
}


void GIF::FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow)
{
	for (int i = 0; i < mNumDOF; i++) {
		int rowInKKT = conesConstraintsStartRow + i;
		int indexOfConeInSystem = mIndexOfHEinSystem[mConesAndMetaMap[i]->halfedge()->index()];
		tripletListValues.push_back(Eigen::Triplet<double>(indexOfConeInSystem, rowInKKT, 1));
	}
}


void GIF::SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow)
{
	//set diagonal elements
	for (int i = 0; i < mSizeOfMatrix; i++)
		tripletListValues.push_back(Eigen::Triplet<double>(i, i, 0));
	for (int i = 0; i != mNonBoundaryVerticesNearConesVerticesVec.size(); i++) {
		int row = mNonBoundaryVerticesNearConesVerticesVec[i];
		int row2 = row + mSizeOfSystemVar;
		for (int j = 0; j < mNumDOF; j++) {
			int col = j + conesConstraintsStartRow;// because we want the cols in the upper right block of the kkt matrix
			tripletListValues.push_back(Eigen::Triplet<double>(row, col, 0));
		}
	}
}


bool GIF::constructHarmonicBasisAndSendToMATLAB(int conesConstraintsStartRow)
{
	//**********calculate harmonic basis in PARDISO************
	bool res;
	res = calculateHarmonicBasisInPARDISO(conesConstraintsStartRow);
	if (!res) return false;
	MatlabInterface::GetEngine().EvalToCout("createJmatrixForOriginalMeshWithScaffold");
	return true;
}

bool GIF::calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow)
{
	//**********calculate harmonic basis in PARDISO************

	int	mtype = -2;//real and symmetric indefinite
	mPardisoSolver = (PardisoLinearSolver(mtype));
	bool res = mPardisoSolver.init(true);
	if (!res) return false;
	res = mPardisoSolver.createPardisoFormatMatrix(mKKtEigen);
	if (!res) return false;
	CGAL::Timer pardiso_timer;
	pardiso_timer.start();
	res = mPardisoSolver.preprocess();
	pardiso_timer.stop();
	std::cout << "Symbolic and Numerical factorization: " << pardiso_timer.time() << "\n";
	pardiso_timer.reset();
	pardiso_timer.start();
	if (!res) return false;
	GMMCompressed1RowMatrix selectiveInvSol(mKKtEigen.rows(), mKKtEigen.cols());
	res = mPardisoSolver.selectiveInverse();
	pardiso_timer.stop();
	std::cout << "Selective Inverse: " << pardiso_timer.time() << "\n";
	if (!res) return false;
	res = mPardisoSolver.getMatrixInGMMformat(selectiveInvSol);
	if (!res) return false;
	//********** Create Harmonic basis matrix for all vertices - **********

	unordered_set<int>::const_iterator itr;
	for (int i = 0; i != mNonBoundaryVerticesNearConesVerticesVec.size(); i++) {
		int row = mNonBoundaryVerticesNearConesVerticesVec[i];
		for (int j = 0; j < mNumDOF; j++) {
			int col = j + conesConstraintsStartRow;

			double realPart = selectiveInvSol(row, col);
			mHarmonicBasisInAllV(row, j) = realPart;
		}
	}
	for (int i = 0; i < mNumDOF; i++) {
		int rowInKKT = i;
		int indexOfConeInSystem = mIndexOfHEinSystem[mConesAndMetaMap[i]->halfedge()->index()];
		mHarmonicBasisInAllV(indexOfConeInSystem, rowInKKT) = 1.0;
	}
	MatlabGMMDataExchange::SetEngineSparseMatrix("HarmonicBasis", mHarmonicBasisInAllV);
	return true;
}


bool GIF::constructScaffoldHarmonicBasisAndSendToMATLAB()
{
	//**********calculate harmonic basis in PARDISO************
	bool res = calculateScaffoldHarmonicBasis();
	if (!res) return false;

	//******* Calculate frames and initialize matrices in MATLAB **********
	MatlabInterface::GetEngine().EvalToCout("createJmatrixScaffold");
	return true;
}

bool GIF::calculateScaffoldHarmonicBasis()
{

	GMMSparseComplexRowMatrix harmonicBasisInAllScaffoldV;
	harmonicBasisInAllScaffoldV.resize(mSizeOfScaffoldSystemVar, mNumScaffoldDOF);
	for (int i = 0; i < mNumScaffoldDOF; i++) {
		int rowInKKT = i;
		int indexOfConeInSystem = mIndexOfHEinScaffoldSystem[mScaffoldConesAndMetaMap[i]->halfedge()->index()];
		harmonicBasisInAllScaffoldV(indexOfConeInSystem, rowInKKT) = 1.0;
	}
	MatlabGMMDataExchange::SetEngineSparseMatrix("ScaffoldHarmonicBasis", harmonicBasisInAllScaffoldV);
	return true;
}


bool GIF::getTutteInitialValue(CGAL_Borders& borders)
{
	GMMDenseColMatrix initialValue(mNumDOF * 2, 1);
	mInitialUVsOfCones.clear();
	mInitialUVsOfCones.resize(mNumDOF);

	int mainBorderIndex = 0;
	Vertex_handle mainBorderVertex = borders.vertex(0, 0);
	int mainBorderVertexIndex = mainBorderVertex->index();

	std::vector<Halfedge_handle> pathMainBorder;
	borders.getBorderHalfEdges(mainBorderVertex, pathMainBorder);
	int min_vertex_index = 0;
	for (int j = 0; j < pathMainBorder.size(); j++)
		if (pathMainBorder[j]->vertex()->index() < pathMainBorder[min_vertex_index]->vertex()->index())
			min_vertex_index = j;
	std::rotate(pathMainBorder.begin(), pathMainBorder.begin() + (min_vertex_index + 1), pathMainBorder.end());
	int numVerticesMainBorder = pathMainBorder.size();

	double totalBorderLength = 0.0;
	for (int i = 0; i < numVerticesMainBorder; i++)
	{
		double edgeLength = pathMainBorder[i]->length();
		totalBorderLength += edgeLength;
	}
	assert(totalBorderLength > 0.0);
	mRadiusOfInitialTutte = 1.0;
	double totalAngle = 0.0;
	double angle = 0.0;
	for (int i = numVerticesMainBorder - 1; i >= 0; --i)
	{
		totalAngle += angle;
		if (mIsMeta[pathMainBorder[i]->vertex()->index()]) {
			double x = cos(totalAngle);
			double y = sin(totalAngle);
			int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[pathMainBorder[i]->vertex()->index()];
			initialValue(indexInUV, 0) = x;
			initialValue(indexInUV + mNumDOF, 0) = y;
			mInitialUVsOfCones[indexInUV] = Complex(x, y);
		}
		double edgeLength = pathMainBorder[i]->length();
		angle = 2.0 * M_PI * edgeLength / totalBorderLength;
	}


	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.UVonCones", initialValue);
	MatlabInterface::GetEngine().EvalToCout("[ GIF.UVonCones, GIF.FixedIndices, GIF.FixedValues ] = fixFirstCone( GIF.UVonCones );");

	return true;
}

bool GIF::runNewtonWithRemeshedScaffold(int conesConstraintsStartRow, GMMDenseColMatrix& RHS)
{
	MatlabInterface::GetEngine().EvalToCout("newtonWithScaffoldPreps;");
	double NewtonSucsess = true;
	GMMDenseColMatrix UVonCones(mNumScaffoldDOF, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones2", UVonCones);
	GMMDenseColMatrix UVonCones1(mNumScaffoldDOF, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones1", UVonCones1);
	int scaffoldConesConstraintsStartRow2;
	int maxRemeshingTimes = 100;
	double curE = 0.0, prevE = 0.0;
	int convex_iteration = 0;
	for (int iter1 = 0; iter1 < maxRemeshingTimes; iter1++)
	{
		MatlabInterface::GetEngine().EvalToCout("iterationOfNewtonWithScaffold;");
		GMMDenseColMatrix UVonCones2(mNumScaffoldDOF, 2);
		MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones1", UVonCones2);
		std::vector<std::complex<double>> internalBorderOfScaffold;
		for (int i = 0; i < mNumDOF; i++)
		{
			internalBorderOfScaffold.push_back(Complex(UVonCones2(mNumScaffoldDOF - i - 1, 0), UVonCones2(mNumScaffoldDOF - i - 1, 1)));
		}
		bool ret = CompilationAccelerator::isSimplePolygon(internalBorderOfScaffold, false);
		if (ret)
		{
			MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones2", UVonCones);
			MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones1", UVonCones1);
		}
		GMMDenseColMatrix arr2(1, 1);
		MatlabGMMDataExchange::GetEngineDenseMatrix("arr1", arr2);

		curE = arr2(0, 0);
		if (((iter1 > 0) && (abs(prevE - curE) / prevE < mOuterRateTerminationCondition)) || (iter1 == maxRemeshingTimes - 1) || !ret || curE <= mEnergyRelatedTerminationCondition)
			break;

		prevE = curE;
		if (iter1 == maxRemeshingTimes - 1)
			break;
		std::vector<std::complex<double>> allScaffoldPoints;
		for (int i = 0; i < mNumScaffoldDOF; i++)
		{
			allScaffoldPoints.push_back(Complex(UVonCones1(i, 0), UVonCones1(i, 1)));
		}

		Mesh curScaffold;
		double outerRadius = 0.0;
		MeshAlgorithms::buildScaffoldWithFixedOuterBoundary(allScaffoldPoints, false, curScaffold, mNumDOF);

		CGAL_Borders borders1(curScaffold);

		Eigen::SparseMatrix<double, Eigen::RowMajor> scaffoldKKtEigen1;

		getScaffoldBordersMapAndSetMetaVertices(borders1, curScaffold.size_of_vertices());
		UpdateIndexOfHEinScaffoldSystem(curScaffold);
		TransferScaffoldBorderFacesAndVerticesToMatlab(curScaffold);


		bool res = constructScaffoldHarmonicBasisAndSendToMATLAB();
		if (!res) {
			cout << "Failed to construct Harmonic Basis for the scaffold" << endl;
			return false;
		}
		int numVerticesScaffold = curScaffold.size_of_vertices();
		GMMDenseColMatrix scaffoldBoundaryVertices(2 * numVerticesScaffold, 1);
		for (int i = 0; i < numVerticesScaffold; i++)
		{
			scaffoldBoundaryVertices(i, 0) = mScaffoldConesAndMetaMap[i]->point().x();
			scaffoldBoundaryVertices(i + numVerticesScaffold, 0) = mScaffoldConesAndMetaMap[i]->point().y();
		}
		MatlabGMMDataExchange::SetEngineDenseMatrix("scaffoldBoundaryVertices", scaffoldBoundaryVertices);
		MatlabInterface::GetEngine().EvalToCout("[ GIF.ScaffoldUVonCones, GIF.FixedIndices, GIF.FixedValues ] = fixFirstCone( scaffoldBoundaryVertices );");


	}
	
	gmm::clear(RHS);
	gmm::resize(RHS, mSizeOfMatrix, 2);
	for (int i = 0; i < mNumDOF; i++)
	{
		int row = conesConstraintsStartRow + i;
		RHS(row, 0) = UVonCones(mNumScaffoldDOF + i, 0);
		RHS(row, 1) = UVonCones(mNumScaffoldDOF + i, 1);
	}

	return (bool)NewtonSucsess;

}


bool GIF::getAllUVsByRHS(GMMDenseColMatrix& RHS)
{
	GMMDenseColMatrix UVs(RHS.nrows(), RHS.ncols());
	bool res1 = mPardisoSolver.solve(RHS, UVs);
	MatlabGMMDataExchange::SetEngineDenseMatrix("array1", UVs);
	MatlabGMMDataExchange::GetEngineDenseMatrix("array1", mUVs);

	return true;
}


//Update the UV's of the halfedges, the input array UVs is ordered such there is a vector for y values and vector for x values ([x , y])
void GIF::updatHalfedgeUVs()
{
	mAllBoundaryUVs.clear();
	mAllBoundaryUVs.resize(mCgalMesh.size_of_border_edges());
	Facet_iterator faceIt = mCgalMesh.facets_begin();
	while (faceIt != mCgalMesh.facets_end())
	{
		Halfedge_handle h = faceIt->halfedge();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1), 0);
		mCgalMesh.vertex(h->vertex()->index())->uvC(Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1)));
		if (h->opposite()->is_border())
		{
			mAllBoundaryUVs[mIndicesOfAllBorderVerticesInUVbyGeneralIndex[h->vertex()->index()]] = Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1));
		}
		h = h->next();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1), 0);
		mCgalMesh.vertex(h->vertex()->index())->uvC(Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1)));
		if (h->opposite()->is_border())
		{
			mAllBoundaryUVs[mIndicesOfAllBorderVerticesInUVbyGeneralIndex[h->vertex()->index()]] = Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1));
		}
		h = h->next();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1), 0);
		if (h->opposite()->is_border())
		{
			mAllBoundaryUVs[mIndicesOfAllBorderVerticesInUVbyGeneralIndex[h->vertex()->index()]] = Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1));
		}
		mCgalMesh.vertex(h->vertex()->index())->uvC(Complex(mUVs(mIndexOfHEinSystem[h->index()], 0), mUVs(mIndexOfHEinSystem[h->index()], 1)));
		faceIt++;
	}
}

void GIF::getResult()
{
	int numFoldovers, numFoldsNearCones, numFoldsNearBorder, numWrongAngles, numWrongConeAngles;
	std::vector<Facet_handle> flippedTriangles;
	checkForFoldovers(flippedTriangles);
	mParameters_matrix(14, 0) = flippedTriangles.size();
	mFlipsCountBeforeFix.clear();
	mFlipsCountBeforeFix.resize(1, 1);
	mFlipsCountAfterFix.clear();
	mFlipsCountAfterFix.resize(1, 1);
	if (!mAllowMVCFix)
	{
		mFlipsCountBeforeFix(0, 0) = flippedTriangles.size();
		cout << mFlipsCountBeforeFix(0, 0) << " flips! not applying fixing method anyway because of user's requirement!" << endl;
		mParameters_matrix(13, 0) = flippedTriangles.size();
		return;
	}

	if (flippedTriangles.size() == 0)
	{
		cout << "No foldovers! No need to apply the foldovers fixing method!" << endl;
		mParameters_matrix(13, 0) = 0;
		mFlipsCountBeforeFix(0, 0) = 0;
		mFlipsCountAfterFix(0, 0) = 0;
	}
	else
	{
		cout << "Applying the foldovers fixing method to fix " << flippedTriangles.size() << " foldovers..." << endl;
		mFlipsCountBeforeFix(0, 0) = flippedTriangles.size();
		fixCotFoldoversWithMeanValue(flippedTriangles);
		checkForFoldovers(flippedTriangles);
		if (flippedTriangles.size() == 0)
			cout << "All foldovers were fixed succesfully!" << endl;
		else
			cout << flippedTriangles.size() << " foldovers remained after the foldovers fixing method was applied" << endl; //this means that the foldovers fixing method didn't work. This is not supposed to happen, because our method has theoretical guarantees
		mParameters_matrix(13, 0) = 1;
		mFlipsCountAfterFix(0, 0) = flippedTriangles.size();
	}

}

void GIF::checkForFoldovers(std::vector<Facet_handle>& flippedTriangles)
{
	flippedTriangles.clear();

	Point_2 a, b, c;
	Facet_iterator faceIt = mCgalMesh.facets_begin();
	while (faceIt != mCgalMesh.facets_end())
	{
		Halfedge_handle h = faceIt->halfedge();
		a = Point_2(h->uv().x(), h->uv().y());
		h = h->next();
		b = Point_2(h->uv().x(), h->uv().y());
		h = h->next();
		c = Point_2(h->uv().x(), h->uv().y());
		Kernel::Triangle_2 t(a, b, c);
		if (t.orientation() <= 0) {
			flippedTriangles.push_back(faceIt);
		}
		faceIt++;
	}
}

void GIF::openReportWindow()
{
	MatlabGMMDataExchange::GetEngineDenseMatrix("array1", mUVs);
	updatHalfedgeUVs();
	std::vector<Facet_handle> flippedTriangles;
	checkForFoldovers(flippedTriangles);
	mParameters_matrix(12, 0) = flippedTriangles.size();

	bool isFoldoversFree = (flippedTriangles.size() > 0);
	bool isSimple = CompilationAccelerator::isSimplePolygon(mAllBoundaryUVs, false);
	if (isFoldoversFree == 0 && isSimple) {
		cout << "The final result is globally injective!" << endl;
		mParameters_matrix(11, 0) = 1;
	}
	else
	{
		mParameters_matrix(11, 0) = 0;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("GIF.parametersMatrix", mParameters_matrix);
	MatlabInterface::GetEngine().EvalToCout("compute_final_energy;");
	MatlabInterface::GetEngine().EvalToCout("GIF_report_window(GIF);");
}

void GIF::fixCotFoldoversWithMeanValue(std::vector<Facet_handle> flippedTriangles)
{
	CGAL::Timer fixingTime1;
	fixingTime1.start();
	GMMDenseColMatrix UVonCones1(mNumScaffoldDOF, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("GIF.UVonCones1", UVonCones1);
	std::vector<std::complex<double>> internalPolygonPoints1;
	std::vector<int> indicesInGluedMeshByIndicesInScaffold;
	double outerRadiusOfScaffold = 0.0;
	for (int i = 0; i < (mNumScaffoldDOF - mNumDOF); i++)
	{
		Complex c1 = Complex(UVonCones1(i, 0), UVonCones1(i, 1));
		indicesInGluedMeshByIndicesInScaffold.push_back(mCgalMesh.size_of_vertices() + i);
	}

	for (int i = 0; i < (mNumScaffoldDOF - mNumDOF); i++)
	{
		Complex c1 = 2.0 * Complex(UVonCones1(i, 0), UVonCones1(i, 1));
		internalPolygonPoints1.push_back(c1);
	}

	for (int i = 0; i < mCgalMesh.size_of_border_edges(); i++)
	{
		internalPolygonPoints1.push_back(mAllBoundaryUVs[i]);
		indicesInGluedMeshByIndicesInScaffold.push_back(mGeneralIndicesUVbyIndicesOfAllBorderVertices[i]);
	}

	Mesh scaffold1;
	double outerRadius = 0.0;

	MeshAlgorithms::buildScaffoldWithFixedOuterBoundary(internalPolygonPoints1, false, scaffold1, mCgalMesh.size_of_border_edges());

	std::vector<Point_3> vertices1;
	std::vector<unsigned int> faceIndices1;
	std::vector<unsigned int> vertexCount1;
	int numVertices1 = mCgalMesh.size_of_vertices();
	int numFaces1 = mCgalMesh.size_of_facets();
	int numFaces2 = scaffold1.size_of_facets();
	GMMDenseColMatrix vertices2(numVertices1, 3);
	GMMDenseColMatrix faces2(numFaces1, 3);
	GMMDenseColMatrix uvsBeforeFixing(numVertices1, 2);
	GMMDenseColMatrix uvsAfterFixing(numVertices1, 2);


	for (int i = 0; i < numVertices1; i++) //inserting the vertices of the original mesh
	{
		Vertex_handle vi = mCgalMesh.vertex(i);
		Complex c1 = vi->uvC();
		Mesh::Point_3 p(c1.real(), c1.imag(), 0.0);
		uvsBeforeFixing(i, 0) = c1.real();
		uvsBeforeFixing(i, 1) = c1.imag();
		vertices1.push_back(p);
	}
	for (int i = 0; i < (mNumScaffoldDOF - mNumDOF); i++)
	{
		Mesh::Point_3 p(scaffold1.vertex(i)->point().x(), scaffold1.vertex(i)->point().y(), 0.0);
		vertices1.push_back(p);
	}
	for (int i1 = 0; i1 < numFaces1; i1++)
	{
		Facet_const_handle f = mCgalMesh.face(i1);

		Halfedge_around_facet_const_circulator h = f->facet_begin();
		const Halfedge_around_facet_const_circulator hEnd = h;

		int numVerticesInCurrentFace = 0;
		CGAL_For_all(h, hEnd)
		{
			faceIndices1.push_back(h->vertex()->index());
			numVerticesInCurrentFace++;
		}
		vertexCount1.push_back(numVerticesInCurrentFace);
	}

	for (int i1 = 0; i1 < numFaces2; i1++)
	{
		Facet_const_handle f = scaffold1.face(i1);

		Halfedge_around_facet_const_circulator h = f->facet_begin();
		const Halfedge_around_facet_const_circulator hEnd = h;

		int numVerticesInCurrentFace = 0;

		CGAL_For_all(h, hEnd)
		{
			faceIndices1.push_back(indicesInGluedMeshByIndicesInScaffold[h->vertex()->index()]);
			numVerticesInCurrentFace++;
		}
		vertexCount1.push_back(numVerticesInCurrentFace);
	}
	Mesh gluedMesh;

	MeshAlgorithms::buildTriangleMesh(vertices1, faceIndices1, gluedMesh);
	std::vector<Eigen::Triplet<double>> tripletListValues1, tripletListValues2, tripletListValues3;
	int gluedMeshBoundaryVerticesCount = gluedMesh.size_of_border_edges();
	int gluedMeshInternalVerticesCount = gluedMesh.size_of_vertices() - gluedMeshBoundaryVerticesCount;
	tripletListValues1.reserve(7 * gluedMesh.size_of_vertices());
	tripletListValues2.reserve(3 * gluedMeshBoundaryVerticesCount);
	//fill the matrix by going on the triangles one by one
	std::vector<int> internalIndices(gluedMesh.size_of_vertices());
	std::vector<int> boundaryIndices(gluedMesh.size_of_vertices());
	GMMDenseColMatrix RHS(gluedMeshBoundaryVerticesCount, 2);
	int num1 = 0, num2 = 0;
	for_each_vertex(vi, gluedMesh)
	{
		if (!vi->is_border())
		{
			internalIndices[vi->index()] = num1;
			num1++;
		}
		if (vi->is_border())
		{
			boundaryIndices[vi->index()] = num2;
			RHS(num2, 0) = -vi->point().x();
			RHS(num2, 1) = -vi->point().y();
			num2++;
		}
	}

	int r = 0; //num free vertices
	int k = 0; //num fixed/dependent vertices
	std::vector<Eigen::Triplet<double>> tripletListValues;
	tripletListValues.reserve(6 * gluedMesh.size_of_vertices() + mNumDOF);
	int n = gluedMesh.countVertexStates(r, k);
	assert(n == r + k);


	std::vector<double> curVertexRow;
	//fill the matrix by going over vertices
	for_each_vertex(v, gluedMesh)
	{
		int i = v->index();
		if (!v->is_border())
		{

			Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
			const Mesh::Halfedge_around_vertex_circulator hEnd = h;
			int numNeighbors = 0;
			CGAL_For_all(h, hEnd)
				numNeighbors++;
			curVertexRow.clear();
			curVertexRow.resize(numNeighbors);
			int neighborOneRingIndex = 0;
			double sumWeights = 0.0;
			CGAL_For_all(h, hEnd)
			{
				double w = h->opposite()->meanValueWeight(false);
				curVertexRow[neighborOneRingIndex] = w;
				neighborOneRingIndex++;
				sumWeights += w;
			}

			neighborOneRingIndex = 0;

			//normalize the weights
			//this is not mandatory since the system should be full ranked, but may improve numerical stability
			CGAL_For_all(h, hEnd)
			{
				int neighborVertexIndex = h->opposite()->vertex()->index();
				double w = curVertexRow[neighborOneRingIndex];
				neighborOneRingIndex++;
				if (!h->opposite()->vertex()->is_border())
					tripletListValues1.push_back(Eigen::Triplet<double>(internalIndices[i], internalIndices[neighborVertexIndex], -w / sumWeights));
				else
					tripletListValues2.push_back(Eigen::Triplet<double>(internalIndices[i], boundaryIndices[neighborVertexIndex], -w / sumWeights));
			}
			tripletListValues1.push_back(Eigen::Triplet<double>(internalIndices[i], internalIndices[i], 1.0));
		}
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> laplacian1;
	laplacian1.resize(gluedMeshInternalVerticesCount, gluedMeshInternalVerticesCount);
	laplacian1.setFromTriplets(tripletListValues1.begin(), tripletListValues1.end());
	laplacian1.makeCompressed();

	Eigen::SparseMatrix<double, Eigen::RowMajor> mat_L12;
	mat_L12.resize(gluedMeshInternalVerticesCount, gluedMeshBoundaryVerticesCount);
	mat_L12.setFromTriplets(tripletListValues2.begin(), tripletListValues2.end());
	mat_L12.makeCompressed();
	const Eigen::SparseMatrix<double> M = mat_L12;
	GMMSparseRowMatrix mat_L12_mean_value(M.rows(), M.cols());

	for (unsigned k = 0; k < (unsigned)M.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
		{
			mat_L12_mean_value(it.row(), it.col()) = it.value();
		}
	}
	GMMDenseColMatrix rhs_elimination(gluedMeshInternalVerticesCount, 2);
	gmm::mult(mat_L12_mean_value, RHS, rhs_elimination);
	GMMDenseColMatrix solution(gluedMeshInternalVerticesCount, 2);
	int	mtype = 11;
	PardisoLinearSolver pardisoSolver(mtype);
	bool res = pardisoSolver.init(false);
	res = pardisoSolver.createPardisoFormatMatrix(laplacian1);
	res = pardisoSolver.preprocess();
	pardisoSolver.solve(rhs_elimination, solution);
	mParameters_matrix(4, 0) = mParameters_matrix(4, 0) + fixingTime1.time();
	num1 = 0;
	for_each_vertex(vi, gluedMesh)
	{
		if (!vi->is_border())
		{
			double x = solution(num1, 0);
			double y = solution(num1, 1);
			vi->uvC(Complex(x, y));
			num1++;
		}
		if (vi->is_border())
		{
			double x = vi->point().x();
			double y = vi->point().y();
			vi->uvC(Complex(x, y));
		}
	}
	for (int i = 0; i < mCgalMesh.size_of_vertices(); i++)
	{
		mCgalMesh.vertex(i)->uvC(gluedMesh.vertex(i)->uvC());
		Complex c1 = gluedMesh.vertex(i)->uvC();
		uvsAfterFixing(i, 0) = c1.real();
		uvsAfterFixing(i, 1) = c1.imag();
	}
	auto heIt = mCgalMesh.halfedges_begin();
	while (heIt != mCgalMesh.halfedges_end())
	{
		heIt->uv() = Point_3(heIt->vertex()->uvC().real(), heIt->vertex()->uvC().imag(), 0);
		mUVs(mIndexOfHEinSystem[heIt->index()], 0) = heIt->vertex()->uvC().real();
		mUVs(mIndexOfHEinSystem[heIt->index()], 1) = heIt->vertex()->uvC().imag();
		heIt++;
	}
	MatlabInterface::GetEngine().EvalToCout("clear array1;");
	MatlabGMMDataExchange::SetEngineDenseMatrix("array1", mUVs);
}
