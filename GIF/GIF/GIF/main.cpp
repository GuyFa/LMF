#include <direct.h>
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "GIF.h"

#define BUFFER_SIZE 4000
void setMatlabPath();

void main(int argc, char* argv[])
{
	if (argc < 2 && argc > 7)
	{
		cout << "You should give 11 - 12 arguments!" << endl;
		return;
	}

	std::string obj_path = argv[1];
	std::string uv_mat_path = argv[2];
	bool run_IDT = strcmp(argv[3], "True") == 0;
	bool allow_MVC_fix = strcmp(argv[4], "True") == 0;
	int interior_faces = stoi(argv[5]);
	int boundary_segment_size = stoi(argv[6]);
	double curvature_meta_vertices_rate = stod(argv[7]);
	double outer_termination_condition_rate = stod(argv[8]);
	double energy_related_termination_condition = stod(argv[9]);
	std::string edge_lengths_mat_path = (argc == 11) ? argv[10] : "";

	setMatlabPath();

	int methodIndex;
	std::string objPath, vfPath;

	GIF program;
	program.run(obj_path, uv_mat_path, run_IDT, allow_MVC_fix, edge_lengths_mat_path, interior_faces, boundary_segment_size, curvature_meta_vertices_rate, outer_termination_condition_rate, energy_related_termination_condition);
	//cin.get();
}

void setMatlabPath()
{
	char path[BUFFER_SIZE];
	_fullpath(path, "..\\", BUFFER_SIZE);
	int index = strlen(path);
	//path[index - 4] = '\0';
	std::string matlabPath(path);
	matlabPath += "MatlabScripts\\";
	MatlabInterface::GetEngine().AddScriptPath(matlabPath.c_str());
}