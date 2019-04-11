// (c) Jim Carstens 2019

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <assert.h>

using std::cout; using std::endl;
using std::ifstream; using std::ofstream;
using std::to_string;
using std::mt19937; using std::uniform_int_distribution; using std::uniform_real_distribution;

const double M_PI = 3.141592653589793238463;
const int NDIM = 3;
const int N = 520;
/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;
double packing_fraction = 0.5;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "fcc_108.dat";
double temperature = 5.;
/* Simulation variables */

double rCut, rCut2, eCut;
const double mass = 1.;
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double v[N][NDIM];
double box[NDIM];

double f[N][NDIM];
double en;

mt19937 _rng;

void prAr(double arr[], int n) {
	cout << '{';
	for (int i = 0; i < n - 1; i++) {
		cout << arr[i];
	}
	cout << arr[n - 1] << endl;
}
void prAr2D(double arr[][NDIM], int n, int dim) {
	cout << '{';
	for (int i = 0; i < n - 1; i++) {
		cout << *arr[i];
		if ((i + 1) % dim == 0) {
			cout << endl;
		}
	}
	cout << arr[n - 1] << endl;
}


void read_data(const char* filename) {
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

void write_data(int step) {
	ofstream myfile;
	myfile.open("coords_step" + to_string(step) + ".output");
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

void set_packing_fraction(void) {
	double volume = 1.0;
	int d, n;
	for (d = 0; d < NDIM; ++d) { volume *= box[d]; }
	double target_volume = (n_particles * particle_volume) / packing_fraction;
	double scale_factor = pow(target_volume / volume, 1.0 / NDIM);
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) { r[n][d] *= scale_factor; }
	}
	for (d = 0; d < NDIM; ++d) { box[d] *= scale_factor; }
}

void init() {
	double vCOM[NDIM] = { 0 };
	double sumV2 = 0.;
	//calculate COM momentum
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] = uniform_real_distribution<double>(-1., 1.)(_rng);
			vCOM[d] += v[i][d] / (double)n_particles;
		}
	}
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			//subtract COM momentum
			v[i][d] -= vCOM[d];		
			//calculate v^2 / N
			sumV2 += (v[i][d] * v[i][d]) / (double)n_particles;
		}
	}
	//calculate fs
	double fs = sqrt(3. * temperature / sumV2);
	//scale vi
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] *= fs;
		}
	}
}

void force() {
	//re-initialise energy and forces
	en = 0.;
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			f[i][d] = 0;
		}
	}
	//calculate forces
	for (int i = 0; i < n_particles - 1; i++) {
		for (int j = i + 1; j < n_particles; j++) {
			double r_ij2 = 0.;
			double r_ij[NDIM] = { 0. };
			for (int d = 0; d < NDIM; d++) {
				double r_ij_d = r[i][d] - r[j][d];
				r_ij_d -= box[d] * (int)(2.0 * r_ij_d / box[d]);//periodic BC
				r_ij2 += r_ij_d;
				r_ij[d] = r_ij_d;
			}
			if (r_ij2 < rCut2) {
				double r2i = 1. / r_ij2;
				double r6i = pow(r2i, 3.);
				double ff = 48. * r2i * r6i * (r6i - 0.5);
				for (int d = 0; d < NDIM; d++) {
					f[i][d] += ff * r_ij[d];
					f[j][d] -= ff * r_ij[d];
					en += 4. * r6i * (r6i - 1.) - eCut;
				}
			}
		}
	}
}


void integrate() {
	double sumV = 0;
	double sumV2 = 0;
	double xm[N][NDIM];
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			double xx = 2. * r[i][d] - xm[i][d] + delta * delta * f[i][d];
			double vi = (xx - xm[i][d]) / (2. * delta);
			sumV2 += vi * vi;
			xm[i][d] = r[i][d];
			r[i][d] -= xx;
			//periodic BC
			r[i][d] -= floor(r[i][d] / box[d]) * box[d];

		}
	}
	temperature = sumV2 / (3. * n_particles);
}

void initCutoff() {
	rCut = 0.5 * box[0];
	rCut2 = rCut * rCut;
	double rCutminus6 = 1 / (pow(rCut, 6.));
	eCut = 4. * (rCutminus6 * (rCutminus6 - 1.));

}

int main() {
	radius = 0.5 * diameter;

	if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
	else if (NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
	else { printf("Number of dimensions NDIM = %d, not supported.", NDIM); return 0; }
	read_data(init_filename);
	set_packing_fraction();
	initCutoff();
	mt19937 _rng(std::chrono::steady_clock::now().time_since_epoch().count());
	init();

	double t = 0.0;
	double deltaT = 0.0001;
	double tMax = 1000.;
	double outputT = 100.;
	for (double t = 0.0; t < tMax; t += deltaT) {
		force();
		integrate();

		if ((int)t % output_steps == 0) {
			cout << "t = " << t << ", E = " << en << ", temperature = " << temperature << endl;
			write_data((int)(t*10000));
		}
	}

	return 0;
}