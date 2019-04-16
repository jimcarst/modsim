#include "molecularDynamics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
using std::cout; using std::endl;
using std::ifstream; using std::ofstream;

molecularDynamics::molecularDynamics()
{
}


molecularDynamics::~molecularDynamics()
{
}


void molecularDynamics::read_data(const char* filename) {
	ifstream ifs(filename);
	if (ifs.fail()) {
		cout << "File doesn't exist" << endl;
		return;
	}
	double dmin, dmax, dia;
	ifs >> n_particles;
	if (n_particles > 10000) { cout << "Too many particles" << endl; return; }
	else if (n_particles <= 0) { cout << "No particles" << endl; return; }
	for (int d = 0; d < NDIM; ++d) {
		ifs >> dmin;
		ifs >> dmax;
		box[d] = fabs(dmax - dmin);
	}
	for (int n = 0; n < n_particles; ++n) {
		for (int d = 0; d < NDIM; ++d) {
			ifs >> r[n][d];
		}
		ifs >> dia;
		if (dia != diameter) { cout << "Wrong diameter" << endl; return; }
	}
	ifs.close();
}

void molecularDynamics::write_data(int step) {
	ofstream myfile;
	myfile.open("coords_step" + std::to_string(step) + ".output");
	int d, n;
	myfile << n_particles << endl;
	for (d = 0; d < NDIM; ++d) {
		myfile << 0.0 << ' ' << box[d] << endl;
	}
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) {
			myfile << r[n][d] << ' ';
		}
		myfile << diameter << endl;
	}
	myfile.close();
}

void molecularDynamics::set_packing_fraction(void) {
	double volume = 1.0;
	int d, n;
	for (d = 0; d < NDIM; ++d) { volume *= box[d]; }
	const double target_volume = (n_particles * particle_volume) / packing_fraction;
	double scale_factor = pow(target_volume / volume, 1.0 / NDIM);
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) { r[n][d] *= scale_factor; }
	}
	for (d = 0; d < NDIM; ++d) { box[d] *= scale_factor; }
}

void molecularDynamics::init() {
	//double *v_com2;
	double *v_com = new double[NDIM] {0};

	//double v_com[NDIM] = { 0 };
	double sum_v2 = 0.;
	//calculate COM momentum
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] = std::uniform_real_distribution<double>(-1., 1.)(rng);
			v_com[d] += v[i][d] / (double)n_particles;
		}
	}
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			//subtract COM momentum
			v[i][d] -= v_com[d];
			//calculate v^2 / N
			sum_v2 += (v[i][d] * v[i][d]) / (double)n_particles;
		}
	}
	//calculate fs
	double fs = sqrt(3. * temperature / sum_v2);
	//scale vi
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] *= fs;
		}
	}
	delete[] v_com;
}

void molecularDynamics::integrate() {
	//double sumV = 0;
	double sum_v2 = 0;

	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			double rTemp_d = 2. * r[i][d] - r_prev_t[i][d] + delta_t * delta_t * f[i][d];
			v[i][d] = (rTemp_d - r_prev_t[i][d]) / (2. * delta_t);
			sum_v2 += (v[i][d] * v[i][d]);// / n_particles;
			r_prev_t[i][d] = r[i][d];
			r[i][d] = rTemp_d;
			//periodic BC
			r[i][d] -= floor(r[i][d] / box[d]) * box[d];

		}
	}
	temperature = sum_v2 / (3. * n_particles);
	double eTot = (en + 0.5 * sum_v2) / n_particles;//unused!

}


void molecularDynamics::force() {
	//re-initialise energy and forces
	en = 0.;
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			f[i][d] = 0.;
		}
	}
	//calculate forces
	for (int i = 0; i < n_particles - 1; i++) {
		for (int j = i + 1; j < n_particles; j++) {

			double r_ij2 = 0.;
			double* r_ij = new double[NDIM] {0.};

			for (int d = 0; d < NDIM; d++) {
				double r_ij_d = r[i][d] - r[j][d];
				r_ij_d -= box[d] * (int)(2.0 * r_ij_d / box[d]);//periodic BC
				r_ij2 += r_ij_d * r_ij_d;
				r_ij[d] = r_ij_d;
			}
			if (r_ij2 < r_cut2) {
				double r6 = 1.0 / (r_ij2 * r_ij2 * r_ij2);
				double ff = 24.0 * r6 * (2.0 * r6 - 1.0);

				for (int d = 0; d < NDIM; d++) {
					f[i][d] += ff * r_ij[d];
					f[j][d] -= ff * r_ij[d];
					en += 4.0 * r6 * (r6 - 1.0) - e_cut;
				}
			}
			delete[] r_ij;
		}
	}

}//blows up!


void molecularDynamics::init_cutoff() {
	r_cut = 0.5 * box[0];
	r_cut2 = r_cut * r_cut;
	double r_cuti6 = 1 / (pow(r_cut, 6.));
	e_cut = 4. * (r_cuti6 * (r_cuti6 - 1.));

}