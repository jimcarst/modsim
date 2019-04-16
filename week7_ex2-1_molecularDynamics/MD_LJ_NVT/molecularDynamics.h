#pragma once
#include <random>
using std::mt19937;

class molecularDynamics
{
public:
	molecularDynamics();
	~molecularDynamics();

	void read_data(const char* filename);
	void write_data(int step);
	void set_packing_fraction(void);
	void init();
	void integrate();
	void force();
	void init_cutoff();

private:
	const int NDIM = 3;
	const int N = 520;

	int n_particles = 0;
	const double diameter = 1.0;
	double packing_fraction = 0.2;
	const double delta_t = 0.0001;
	double temperature = 100.;
	double radius;
	double r_cut, r_cut2, e_cut;
	double particle_volume;
	double** r;

	double** v;
	double* box;
	double** r_prev_t;
	double** f;
	double en;
	mt19937 rng;

};

