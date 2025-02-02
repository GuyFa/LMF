#pragma once

#include "Energy.h"
#include "EigenTypes.h"

#include <atomic>
#include <functional>
#include <shared_mutex>
using namespace std;

class Solver
{
public:
	Solver();

	int run();
	void stop();
	void get_mesh(MatX2& X);
	void init(const MatX3& V, const MatX3i& F, float& time_Tutte_init);
	void init_from_input(const MatX3& V, const MatX3i& F, const MatX3& init_uv);

	void computeTutte(const MatX3& V, const MatX3i& F, MatX2 &uv_init);
	// Pointer to the energy class
	shared_ptr<Energy> energy;

	// Activity flags
	atomic_bool is_running = {false};
	atomic_bool progressed = {false};

	// Answer from explicit solver after step done
	int ret;

	// Synchronization functions used by the wrapper
	void wait_for_param_slot();
	void release_param_slot();

	// Soup vertices and face indices
	MatX2 uv;
	MatX3i F;

	// External (interface) and internal working mesh
	Vec ext_x, m_x;

	int num_steps = 2147483647;

	bool full_init_done = false;

	// Current energy, gradient and hessian
	double f;
	Vec g;
	SpMat h;


protected:
	// Give the wrapper a chance to intersect gracefully
	void give_param_slot();
	// Updating the mesh after a step has been done
	void update_external_mesh();

	// Descent direction evaluated in step
	Vec p;

	// Function pointers to the full and value-only energy evaluation
	function<void(const Vec&, double&)> eval_f;
	function<void(const Vec&, double&, Vec&, SpMat&)> eval_fgh;
	
	
	// Synchronization structures
	atomic_bool params_ready_to_update = {false};
	atomic_bool wait_for_param_update = {false};
	atomic_bool a_parameter_was_updated = {false};
	atomic_bool halt = {false};
	
	// Mutex needed in lbfgs - thus protected
	unique_ptr<shared_timed_mutex> mesh_mutex;

	// pardiso variables
	vector<int> IId, JJd, IIp, JJp, II, JJ;
	vector<double> SSd, SSp, SS;

private:
	virtual int step() = 0;
	virtual double linesearch() = 0;
	virtual bool test_progress() = 0;
	virtual void internal_init() = 0;
	virtual void internal_update_external_mesh() = 0;

	// Mutex stuff
	unique_ptr<mutex> param_mutex;
	unique_ptr<condition_variable> param_cv;
	
};