// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <stdio.h>
#include "Lattice.h"
#include <vector>

std::mt19937 rng;
std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
double temperature, beta;
int J;
const int N = 50;
const int lattice_size = N * N;
//int lattice[N * N] = { 0 };
double prob[18];
double p_add;
Lattice lat(N, N);
int M_total, E_total;

template <class T>
double mean(std::vector<T>& buffer) {
	double sum = 0;
	for (int i = 0; i < buffer.size(); i++) {
		sum += buffer.at(i);
	}
	return sum / buffer.size();
}
template <class T>
double mean_squared(std::vector<T>& buffer) {
	double sum = 0;
	for (int i = 0; i < buffer.size(); i++) {
		sum += buffer.at(i) * buffer.at(i);
	}
	return sum / buffer.size();
}

template <class T>
double variance(std::vector<T>& buffer) {
	double sd = 0;
	int n = buffer.size();
	double mu = mean(buffer);
	for (int i = 0; i < n; i++) {
		sd += (buffer.at(i) - mu) * (buffer.at(i) - mu);
	}
	sd /= (n - 1.);
	return sd;
}


void print() {
	for (int i = 0; i < N; i++) cout << "=";
	cout << "\n";
	for (int i = 0; i < lattice_size; i++) {
		//cout << std::showpos << lat.at(i) << ",";
		cout << (lat.at(i) == 1 ? (char)250 : (char)254);
		if ((i + 1) % N == 0) {
			cout << "\n";
		}
	}
	for (int i = 0; i < N; i++) cout << "=";
	cout << "\n";
}

void write_data(int step) {
	char buffer[128];
	snprintf(buffer, 128, "configurations/conf_step%d.output", step);
	std::ofstream file;
	file.open(buffer);
	for (int x = 0; x < N; x++) {
		for (int y = 0; y < N - 1; y++) {
			file << lat(x, y) << ' ';
		}
		file << lat(x, N - 1) << "\n";
	}
	file.close();
}
void read_data(const char* filename) {
	std::ifstream ifs(filename);
	if (ifs.fail()) {
		cout << "File doesn't exist" << endl;
		return;
	}
	int value, count(0);
	while (ifs >> value) {
		cout << value << ", ";
		lat.set_at(count, value);
		count++;
	}
}

int magnetisation() {
	int M = 0;
	for (int i = 0; i < lattice_size; i++) {
		M += lat.at(i);
	}
	return M;//per latice
}
int energy() {
	int E = 0;
	for (int x = 0; x < N; x++) {
		for (int y = 0; y < N; y++) {
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {
					if (!(dx == 0 && dy == 0)) {
						int neighbour_x = x + dx;
						int neighbour_y = y + dy;
						if (neighbour_x == N) {
							neighbour_x = 0;
						}
						if (neighbour_x == -1) {
							neighbour_x = N - 1;
						}
						if (neighbour_y == N) {
							neighbour_y = 0;
						}
						if (neighbour_y == -1) {
							neighbour_y = N - 1;
						}
						E += lat(x, y) * lat(neighbour_x, neighbour_y);
					}
				}
			}
		}
	}
	E *= -J;
	return E;
}

void set_probabilities() {
	for (int i = 2; i < 18; i += 2) {
		prob[i] = exp(-beta * i);
	}
	p_add = 1 - exp(-2.0 * beta * J);
}
void init_infinite_temperature() {
	for (int x = 0; x < N; x++) {
		for (int y = 0; y < N; y++) {
			if (random_zero_one(rng) <= 0.5) {
				lat(x, y) = -1;
			}
			else {
				lat(x, y) = 1;
			}
		}
	}
}

void flip(int x, int y) {
	lat.set(x, y, -lat(x, y));
}

void sweep() {
	for (int i = 0; i < lattice_size; i++) {
		int random_site = std::uniform_int_distribution<int>(0, lattice_size - 1)(rng);
		int random_x = random_site % N;
		int random_y = random_site / N;
		int sum = 0;
		//for (int n : {i - N - 1, i - N, i - N + 1, i - 1, i + 1, i + N - 1, i + N, i + N + 1}) {}
		//above could be faster, but checking periodic boundary cond's becomes messy.
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (dx == 0 && dy == 0) continue;
				int neighbour_x = random_x + dx;
				int neighbour_y = random_y + dy;
				if (neighbour_x == N) {
					neighbour_x = 0;
				}
				if (neighbour_x == -1) {
					neighbour_x = N - 1;
				}
				if (neighbour_y == N) {
					neighbour_y = 0;
				}
				if (neighbour_y == -1) {
					neighbour_y = N - 1;
				}
				sum += lat(neighbour_x, neighbour_y);
			}
		}
		int delta_E = 2 * J * lat(random_x, random_y) * sum;
		if (delta_E <= 0) {
			flip(random_x, random_y);
			E_total += 2 * delta_E;
			M_total += 2 * lat(random_x, random_y);
		}
		else if (random_zero_one(rng) < prob[delta_E]) {
			flip(random_x, random_y);
			E_total += 2 * delta_E;
			M_total += 2 * lat(random_x, random_y);
		}
	}
}
void wolff_step() {
	int stack[lattice_size];
	int seed = std::uniform_int_distribution<int>(0, lattice_size - 1)(rng);
	int seed_x = seed % N;
	int seed_y = seed / N;
	stack[0] = seed;
	int sp = 1; // sp is the size of the LIFO buffer
	int oldspin = lat.at(seed);
	int newspin = -lat.at(seed);
	flip(seed_x, seed_y);
	while (sp) {
		int current = stack[--sp];
		int current_x = current % N;
		int current_y = current / N;
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (dx == 0 && dy == 0) continue;
				int neighbour_x = current_x + dx;
				int neighbour_y = current_y + dy;
				//periodic BC's
				if (neighbour_x == N) neighbour_x = 0;
				if (neighbour_x == -1) neighbour_x = N - 1;
				if (neighbour_y == N) neighbour_y = 0;
				if (neighbour_y == -1) neighbour_y = N - 1;
				if (lat(neighbour_x, neighbour_y) == oldspin) {
					//p_add = 1-exp(-2*beta*J)
					if (random_zero_one(rng) < p_add) {
						//add to cluster
						stack[sp++] = neighbour_x + N * neighbour_y;
						flip(neighbour_x, neighbour_y);
					}
				}
			}
		}
	}
}

int main() {
	const int mc_steps = 1100000;
	const int output_steps = 1000;
	//lat.init(1);
	init_infinite_temperature();
	temperature = 5.0;
	beta = 1. / temperature;
	set_probabilities();
	J = 1;	
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	M_total = magnetisation();
	E_total = energy();
	std::ofstream measure_file;
	measure_file.open("measurements_M_E.output");
	std::ofstream plot_file;
	plot_file.open("plotdata.output");
	measure_file << M_total << ',' << E_total << "\n";
	int measure_series = 0;
	int measure_time = 0;
	std::vector<double> e_measurement;
	std::vector<double> m_measurement;
	for (int i = 0; i < mc_steps; ++i) {
		//sweep();
		wolff_step();
		if (i % 10 == 0) write_data(i);
		measure_time++;
		measure_file << M_total << ',' << E_total << "\n";
		/*if (measure_time > 1000) {
			measure_file.close();
			plot_file.close();
			return 0;
		}*/

		if (measure_time > 1000 && measure_time % 100 == 0) {
			e_measurement.push_back((double)E_total / lattice_size);
			m_measurement.push_back((double)M_total / lattice_size);
		}

		if (measure_time > 11000) {//1000+100*100
			double e_mean = mean(e_measurement);
			double m_mean = mean(m_measurement);
			//cout << "MEAN E = " << e_mean << " variance = " << variance(e_measurement);
			//cout << ", mean M = " << m_mean << ", variance = " << variance(m_measurement) << endl;
			double heat_cap = beta * beta * (mean_squared(e_measurement) - e_mean * e_mean) * (double)lattice_size;
			double magnetic_sus = beta * (double)lattice_size * (mean_squared(m_measurement) - m_mean * m_mean);
			//cout << "hc = " << heat_cap << ", msus = " << magnetic_sus << endl;
			plot_file << temperature << ',' << m_mean << ',' << variance(m_measurement) << ',' << e_mean;
			plot_file << ',' << variance(e_measurement) << ',' << heat_cap << ',' << magnetic_sus << "\n";

			cout << "WRITE (T, m, e, heat_cap, mag_sus) " << temperature << ", " << m_mean << ", " << e_mean;
			cout << ", " << heat_cap << ", " << magnetic_sus << "\n";
			measure_time = 0;
			e_measurement.clear();
			m_measurement.clear();
			temperature -= 0.1;
			beta = 1. / temperature;
			set_probabilities();
		}

		//if (i % output_steps == 0) {			
			//cout << "m = " << (double)M_total/lattice_size << ", E = " << (double)E_total/lattice_size << ", \n";
		//}
		//if (i % 20000 == 0) {
		//	temperature -= 0.3;
		//	beta = 1. / temperature;
		//	set_probabilities();
		//	cout << "T = " << temperature << ", " << (char)225 << " = " << beta << ", m = ";
		//	cout << (double)M_total / lattice_size << ", E/N = " << (double)E_total / lattice_size;
		//	cout << ", E = " << E_total << ", E_real = " << energy() << "\n";
		//}
	}
	write_data(1);

	measure_file.close();
	plot_file.close();
	return 0;
}