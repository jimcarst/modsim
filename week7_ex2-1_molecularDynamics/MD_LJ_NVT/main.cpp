// (c) Jim Carstens 2019

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
using std::cout; using std::endl;

const double M_PI = 3.141592653589793238463;
const int NDIM = 3;
const int N = 520;

//Step_variables
const double diameter = 1.0;
double temperature = 10.;

//MD_Variables
double r_cut, r_cut2, e_cut;
const double delta_t = 0.00001;

int n_particles = 0;
double particle_volume;
double e_potential, e_total;

double r[N][NDIM];
double v[N][NDIM];
double box[NDIM];
double r_prev_t[N][NDIM];
double f[N][NDIM];

std::mt19937 rng;

void read_data(const char* filename) {
	std::ifstream ifs(filename);
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
	std::ofstream myfile;
	myfile.open("coords_step" + std::to_string(step) + ".output");
	int d, n;
	myfile << n_particles << "\n";
	for (d = 0; d < NDIM; ++d) {
		myfile << 0.0 << ' ' << box[d] << "\n";
	}
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) {
			myfile << r[n][d] << ' ';
		}
		myfile << diameter << "\n";
	}
	myfile.close();
}

void set_packing_fraction(double packingFraction) {
	double volume = 1.0;
	int d, n;
	for (d = 0; d < NDIM; ++d) { volume *= box[d]; }
	const double target_volume = (n_particles * particle_volume) / packingFraction;
	double scale_factor = pow(target_volume / volume, 1.0 / NDIM);
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) { r[n][d] *= scale_factor; }
	}
	for (d = 0; d < NDIM; ++d) { box[d] *= scale_factor; }
}

void init() {
	double v_com[NDIM] = { 0. };
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
}

void init_cutoff() {
	r_cut = 0.5 * box[0];
	r_cut2 = r_cut * r_cut;
	double r_cuti6 = 1. / (pow(r_cut, 6.));
	e_cut = 4. * (r_cuti6 * (r_cuti6 - 1.));
}

void force() {
	//re-initialise energy and forces
	e_potential = 0.;
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			f[i][d] = 0.;
		}
	}
	//calculate forces
	for (int i = 0; i < n_particles - 1; i++) {
		for (int j = i + 1; j < n_particles; j++) {			
			double r_ij2 = 0.;
			double r_ij[NDIM] = { 0. };
			for (int d = 0; d < NDIM; d++) {
				r_ij[d] = r[i][d] - r[j][d];
				//r_ij[d] -= box[d] * (int)(2.0 * r_ij[d] / box[d]);//Nearest Image Convention
				if (r_ij[d] > 0.5 * box[d]) {
					r_ij[d] -= box[d];
				}
				if (r_ij[d] < -0.5 * box[d]) {
					r_ij[d] += box[d];
				}
				
				r_ij2 += r_ij[d] * r_ij[d];
				//r_ij[d] = r_ij_d;
			}	
			if (r_ij2 < r_cut2) {
				double r6 = 1.0 / (r_ij2 * r_ij2 * r_ij2);
				double ff = 24.0 * r6 * (2.0 * r6 - 1.0);

				for (int d = 0; d < NDIM; d++) {
					f[i][d] += ff * r_ij[d];
					f[j][d] -= ff * r_ij[d];
					e_potential += 1.0 * r6 * (r6 - 1.0) - 0.25*e_cut;
				}
			}
		}
	}
}


void integrate() {
	double sum_v2 = 0.;

	for (int i = 0; i < n_particles; i++) {
		double r_temp[NDIM] = { 0. };
		for (int d = 0; d < NDIM; ++d) {
			r_temp[d] = 2. * r[i][d] - r_prev_t[i][d] + delta_t * delta_t * f[i][d];
		}
		for (int d = 0; d < NDIM; ++d) {			
			v[i][d] = (r_temp[d] - r_prev_t[i][d]) / (2. * delta_t);
			sum_v2 += (v[i][d] * v[i][d]);// / n_particles;
		}
		for (int d = 0; d < NDIM; ++d) {
			r_prev_t[i][d] = r[i][d];
			r[i][d] = r_temp[d];			
			r[i][d] -= floor(r[i][d] / box[d]) * box[d];//periodic BC
		}
	}
	temperature = sum_v2 / (3. * n_particles);
	e_total = (e_potential + 0.5 * sum_v2); //e_total = e_potential + e_kin
	//cout <<"{"<< e_potential << ", " << sum_v2 << "},\n";
}


int main() {
	const int mc_steps = 1000000;
	const int output_steps = 100;
	double packing_fraction = 0.1;
	const char* init_filename = "fcc_108.dat";

	std::ios::sync_with_stdio(false);

	//particle volume
	if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
	else if (NDIM == 2) particle_volume = 0.25 * M_PI * diameter * diameter;
	else { cout << "Number of dimensions NDIM = " << NDIM << ", not supported."; return 0; }

	//set start parameters
	read_data(init_filename);
	set_packing_fraction(packing_fraction);
	init_cutoff();
	init();
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

	double t = 0.0;	

	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			r_prev_t[i][d] = r[i][d]-v[i][d]*delta_t;
		}
	}
	std::ofstream myfile;
	myfile.open("energy.output");

	for (int i = 0; i < mc_steps; i++) {
		t = i * delta_t;
		force();
		integrate();
		if (i % output_steps == 0) {
			cout << "t = " << t << ", E = " << e_total << ", temperature = " << temperature << "\n";
			myfile << t << ',' << e_total << ',' << temperature << "\n";
			write_data((int)(t*10000));
		}
	}
	myfile.close();
	return 0;
}