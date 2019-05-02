// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <iomanip>
using std::cout; using std::endl; using std::setw;
const double M_PI = 3.141592653589793238463;
const int NDIM = 3;
const int N = 520;
//Step_variables
const double diameter = 1.0;
const double temperature = 1.;
double temp_instantaneous;
//MD_Variables
double r_cut, r_cut2, e_cut;
const double delta_t = 0.000001;
double gamma;
int n_particles = 0;
double particle_volume;
double e_potential, e_total, e_kinetic;
double r[N][NDIM];
double r0[N][NDIM];
double v[N][NDIM];
double v0[N][NDIM];
double v_v_autocorrelation;
double mean_square_displacement;
double box[NDIM];
double r_prev_t[N][NDIM];
double f[N][NDIM];
int version;
std::mt19937 rng;
std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
std::uniform_real_distribution<double> random_minusOne_one(-1.0, 1.0);
std::normal_distribution<double> gaussian(0.0, sqrt(temperature));

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
	std::ofstream energy_file;
	energy_file.open("coords_step" + std::to_string(step) + ".output");
	int d, n;
	energy_file << n_particles << "\n";
	for (d = 0; d < NDIM; ++d) {
		energy_file << 0.0 << ' ' << box[d] << "\n";
	}
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) {
			energy_file << r[n][d] << ' ';
		}
		energy_file << diameter << "\n";
	}
	energy_file.close();
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

void init_velocities_positions() {
	double v_com[NDIM] = { 0. };
	double sum_v2 = 0.;
	//calculate COM momentum
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; ++d) {
			v[i][d] = random_minusOne_one(rng);
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
			r_prev_t[i][d] = r[i][d] - v[i][d] * delta_t;
			v0[i][d] = v[i][d];
			r0[i][d] = r[i][d];
		}
	}
}

void init_cutoff() {
	if (version == 0) {//Lennard-Jones
		r_cut = 2.5 * diameter;
	}
	else if (version == 1) {//Weeks Chandler Anderson
		r_cut = pow(2.0, 1.0 / 6.0) * diameter;
	}
		r_cut2 = r_cut * r_cut;
		double r_cuti6 = 1. / (pow(r_cut, 6.));//2.5 sigma
		e_cut = 4. * (r_cuti6 * (r_cuti6 - 1.));
}

void force() {
	//re-initialise energy and forces
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
				//Nearest Image Convention
				if (r_ij[d] > 0.5 * box[d]) {
					r_ij[d] -= box[d];
				}
				if (r_ij[d] < -0.5 * box[d]) {
					r_ij[d] += box[d];
				}
				r_ij2 += r_ij[d] * r_ij[d];
			}
			if (r_ij2 < r_cut2) {
				double r2i = 1.0 / r_ij2;
				double r6i = r2i * r2i * r2i;
				double ff = 48.0 * r2i * r6i * (r6i - 0.5);
				for (int d = 0; d < NDIM; d++) {
					f[i][d] += ff * r_ij[d];
					f[j][d] -= ff * r_ij[d];
				}
			}
		}
	}
}

void integrate() {
	for (int i = 0; i < n_particles; i++) {
		double r_temp[NDIM] = { 0. };
		for (int d = 0; d < NDIM; d++) {
			double random = gaussian(rng);
			r_temp[d] = r[i][d] + sqrt(2. * temperature * delta_t / gamma) * random + f[i][d] * delta_t / gamma;
			r_temp[d] -= floor(r_temp[d] / box[d]) * box[d];//periodic BC
			v[i][d] = sqrt(2. * temperature / gamma) * random + f[i][d] / gamma;
		}
		for (int d = 0; d < NDIM; ++d) {
			double delta_r = r_temp[d] - r[i][d];
			if (delta_r > 0.5 * box[d]) {// Nearest Image Convention
				delta_r -= box[d];
			}
			if (delta_r < -0.5 * box[d]) {
				delta_r += box[d];
			}
		//	v[i][d] = delta_r / delta_t;
			r[i][d] = r_temp[d];
		}
	}
}

void sample() {
	//potential energy
	e_potential = 0.;
	//calculate forces
	double r_ij2;
	for (int i = 0; i < n_particles - 1; i++) {
		for (int j = i + 1; j < n_particles; j++) {
			//double r_ij2 = 0.;
			r_ij2 = 0.;
			double r_ij[NDIM] = { 0. };
			for (int d = 0; d < NDIM; d++) {
				r_ij[d] = r[i][d] - r[j][d];
				if (r_ij[d] > 0.5 * box[d]) {
					r_ij[d] -= box[d];
				}
				if (r_ij[d] < -0.5 * box[d]) {
					r_ij[d] += box[d];
				}
				r_ij2 += r_ij[d] * r_ij[d];
			}
			if (r_ij2 < r_cut2) {
				double r2i = 1.0 / r_ij2;
				double r6i = r2i * r2i * r2i;
				if (version == 0) {
					e_potential += 4.0 * r6i * (r6i - 1.0) - e_cut;
				}
				else if (version == 1) {
					e_potential += 4.0 * r6i * (r6i - 1.0) + 1.;
				}
			}
		}
	}
	//kinetic energy
	double sum_v2 = 0.;
	v_v_autocorrelation = 0.;
	mean_square_displacement = 0.;
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; d++) {
			sum_v2 += v[i][d] * v[i][d];
			v_v_autocorrelation += v0[i][d] * v[i][d] / n_particles;
			mean_square_displacement += (r[i][d] - r0[i][d]) * (r[i][d] - r0[i][d]) / n_particles;
		}
	}
	//temperature
	temp_instantaneous = sum_v2 /(3. * n_particles);//s = 3.

	e_kinetic = 0.5 * sum_v2;
	e_total = (e_potential + e_kinetic); //e_total = e_potential + e_kin
}

int main() {
	const int mc_steps = 2000000;
	const int output_steps = 1000;
	double packing_fraction = 1.0;
	version = 0;//Lennard Jones
	//version = 1;//Weeks Chandler Anderson
	gamma = 6.0 * M_PI * packing_fraction * diameter;
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
	init_velocities_positions();
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

	double t = 0.0;

	std::ofstream energy_file;// , correlation_file;
	energy_file.open("energy.output");
	//correlation_file.open("correlation.output");
	force();
	v_v_autocorrelation = 0.0;
	mean_square_displacement = 0.0;
	double D = 0.0;
	for (int i = 0; i < mc_steps; i++) {
		t = i * delta_t;
		force();
		integrate();
		D += (1. / (double)NDIM) * delta_t * v_v_autocorrelation;
		if (i % output_steps == 0) {
		sample();
			cout << "t = " << t << ", E = " << e_total << ", temp = " << temperature
				<< ", temp_inst = " << temp_instantaneous << ", e_kin = " << e_kinetic
				<< ", e_pot = " << e_potential << "\n";
			//cout << "{" << e_potential << ", " << e_kinetic << "},\n";
			energy_file << t << ',' << e_total << ',' << temp_instantaneous << ',' << e_potential << ',' << e_kinetic << ',' << v_v_autocorrelation << ',' << mean_square_displacement << "\n";
			//correlation_file << v_v_autocorrelation << ',' << mean_square_displacement << "\n";
			write_data((int)(t * 10000));
		}
	}
	cout << "D = " << D << "\n";
	energy_file.close();
	return 0;
}