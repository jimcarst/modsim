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
const int output_steps = 1000;
double packing_fraction = 0.5;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "fcc_108.dat";
double temperature = 10;
/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double v[N][NDIM];
double box[NDIM];


double f[N];
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
	myfile.open("coords_step" + to_string(step) + ".output");//add step
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


int move_particle(void) {
	int rpid = uniform_int_distribution<int>(0, n_particles)(_rng);
	if (rpid >= n_particles || rpid < 0) {//here???
		return 0;
	}
	double new_pos[NDIM];
	int d;
	for (d = 0; d < NDIM; ++d) {
		new_pos[d] = r[rpid][d] + uniform_real_distribution<double>(-delta, delta)(_rng);
		new_pos[d] -= floor(new_pos[d] / box[d]) * box[d];
		//assert(new_pos[d] < box[d]);
		//assert(new_pos[d] >= 0.0);
	}
	for (int n = 0; n < n_particles; ++n) {
		if (n == rpid) continue;
		double dist = 0.0;
		for (d = 0; d < NDIM; ++d) {
			double min_d = new_pos[d] - r[n][d];
			if (min_d > 0.5 * box[d]) {
				min_d -= box[d];
			}
			if (min_d < -0.5 * box[d]) {
				min_d += box[d];
			}
			// assert(min_d <= 0.5 * box[d]);
			// assert(min_d >= -0.5 * box[d]);
			dist += min_d * min_d;
		}
		if (dist <= diameter * diameter) {
			return 0;
		}
	}

	//Accept the move if reached here
	for (d = 0; d < NDIM; ++d) r[rpid][d] = new_pos[d];

	return 1;
}

void set_packing_fraction(void) {
	double volume = 1.0;
	int d, n;
	for (d = 0; d < NDIM; ++d) volume *= box[d];

	double target_volume = (n_particles * particle_volume) / packing_fraction;
	double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
	}
	for (d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

void init() {
	double sumV = 0.;
	double sumV2 = 0.;
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] = uniform_real_distribution<double>(-1., 1.)(_rng);
			sumV += v[i][d];
			sumV2 += v[i][d] * v[i][d];
		}
	}
	sumV /= n_particles;
	sumV2 /= n_particles;
	double fs = sqrt(3 * temperature / sumV2);
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] = (v[i][d] - sumV) * fs;
		}
	}
}

void force() {
	en = 0;
	double eCut = 10;//formula
	double rc = 1;
	double rc2 = rc * rc;
	for (int i = 0; i < n_particles; i++) {
		f[i] = 0;
	}
	for (int i = 0; i < n_particles - 1; i++) {
		for (int j = i + 1; j < n_particles; j++) {
			for (int d = 0; d < NDIM; ++d) {
				double xR = r[i][d] - r[j][d];
				xR -= box[d] * (int)(xR / box[d]);//periodic BC
				double r2 = xR * xR;
				if (r2 < rc2) {
					double r2i = 1 / r2;
					double r6i = pow(r2i, 3.);
					double ff = 48 * r2i * r6i * (r6i - 0.5);
					f[i] += ff * xR;
					f[j] -= ff * xR;
					en += 4 * r6i * (r6i - 1) - eCut;
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
			double xx = 2 * r[i][d] - xm[i][d] + delta * delta * f[i];
			double vi = (xx - xm[i][d]) / (2 * delta);
			sumV2 += vi * vi;
			xm[i][d] = r[i][d];
			r[i][d] -= xx;
			//periodic BC
			r[i][d] -= floor(r[i][d] / box[d]) * box[d];

		}
	}
	temperature = sumV2 / (3. * n_particles);
}

int main() {
	radius = 0.5 * diameter;

	if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
	else if (NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
	else { printf("Number of dimensions NDIM = %d, not supported.", NDIM); return 0; }
	read_data(init_filename);

	set_packing_fraction();
	mt19937 _rng(std::chrono::steady_clock::now().time_since_epoch().count());
	init();

	int accepted = 0;
	double t = 0.0;
	double deltaT = 0.1;
	double tMax = 10;

	for (double t = 0.0; t < tMax; t += deltaT) {
		force();
		integrate();
		cout << "t = " << t << " , e = " << en << endl;
		write_data(t);
	}


	/*
	for (int step = 0; step < mc_steps; ++step) {
		for (int n = 0; n < n_particles; ++n) {
			accepted += move_particle();
		}
		if (step % output_steps == 0) {
			cout << "Step " << step << " Move acceptance: " << (double)accepted / ((double)n_particles * output_steps) << endl;
			accepted = 0;
			write_data(step);
		}*/

	return 0;
}