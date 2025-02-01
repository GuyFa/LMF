#include "Newton.h"
#include <igl\read_triangle_mesh.h>
#include <igl\writeOBJ.h>
#include <time.h>



int main(int argc, char** argv)
{
	unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	omp_set_num_threads(num_threads);

	if (argc < 3)
	{
		cout << "Syntax: Parameterization_cmd.exe <Input file name> <Output file name> <Meta data .txt file name> <target energy> <target time> <Initialization file name>" << endl;
		return false;
	}
	// Tutte is performed anyway, even if target time is smaller
	MatX3 V;
	MatX3i F;
	cout << "Started loading mesh..." << endl;
	if (!igl::read_triangle_mesh(argv[1], V, F))
	{
		cerr << "Failed to load mesh: " << argv[1] << endl;
		return false;
	}

	MatX3 init_V;
	MatX3i init_F;

	if (argc == 7)
		igl::read_triangle_mesh(argv[6], init_V, init_F);

	float target_Esd = 2.0;
	if (argc >= 5)
		target_Esd = stof(argv[4]);

	float target_time = 100000.0;
	if (argc >= 6)
		target_time = stof(argv[5]);

	cout << "Finished loading mesh..." << endl;

	float time_Tutte_init = 0.0;
	long start_init_time = clock();

	Newton solver;
	if (argc <= 6)
		solver.init(V, F, time_Tutte_init);
	else
		solver.init_from_input(V, F, init_V);

	long end_init_time = clock();

	long cur_time = clock();

	float prevF, curF;
	Vec prevX = solver.m_x;

	int fcounter = 0, xcounter = 0, num_conv_iters = 3;
	float ftol = 1e-2, xtol = 1e-10;
	int iter = 0;

	prevF = solver.f / solver.energy->symDirichlet->Area.sum();

	int MAX_ITERS = 1000;
	int iters_count = 0;

	for (iter = 0; iter < MAX_ITERS; iter++)
	{
		solver.step();
		curF = solver.linesearch();

		cur_time = clock();
		if ((cur_time - start_init_time) / 1000.0 >= target_time)
			break;

		iters_count++;


		if (abs(curF - prevF) < ftol * (solver.f + 1))
		{
			if (fcounter >= num_conv_iters)
			{
				cout << "Converged: Change in energy < tol" << endl;
				break;
			}
			else
				fcounter += 1;
		}
		else
			fcounter = 0;

		if ((solver.m_x - prevX).norm() < xtol * (solver.m_x.norm() + 1))
		{
			if (xcounter >= num_conv_iters)
			{
				cout << "Converged: Change in X < tol" << endl;
				break;
			}
			else
				xcounter += 1;
		}
		else
			xcounter = 0;

		if (curF <= target_Esd)
		{
			cout << "Converged: reached target Esd" << endl;
			break;
		}

		Vec tmp = solver.m_x;
		prevX = tmp;
		prevF = solver.f / solver.energy->symDirichlet->Area.sum();

	}

	long end_Newton_time = clock();


	MatX3 CN;
	MatX3i FN;
	solver.uv = Eigen::Map<MatX2>(prevX.data(), prevX.size() / 2, 2);
	MatX3 UVs(V.rows(), 3);
	UVs.leftCols(2) = solver.uv;
	UVs.rightCols(1).setZero();
	if (strcmp(argv[2], "") != 0)
		igl::writeOBJ(argv[2], UVs, F, CN, FN, CN, FN);

	std::string file_path = argv[3];

	// Create an ofstream object
	std::ofstream file(file_path);

	// Check if the file is open
	if (file.is_open()) {
		file << "Total time consumption: " << (end_Newton_time - start_init_time) / 1000.0 << endl;
		file << "Newton iteration count: " << iters_count << std::endl;
		file << "Final E_SD: " << curF << endl;
		file << "Tutte initialization time: " << time_Tutte_init << endl;
		file << "Tutte + Hessian initialization time: " << (end_init_time - start_init_time) / 1000.0 << endl;
		file << "Did Tutte was finished before target time: " << (time_Tutte_init <= target_time) << endl;


		// Close the file
		file.close();

		std::cout << "Text successfully saved to " << file_path << std::endl;
	}
	else {
		// Handle the error if the file couldn't be opened
		std::cerr << "Failed to open file " << file_path << std::endl;
	}
	return 0;
}