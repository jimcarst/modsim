// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <stdio.h>
#include "Lattice.h"


std::mt19937 rng;
std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
double temperature, beta;
int J;
const int N = 50;
const int lattice_size = N * N;
//int lattice[N * N] = { 0 };
double prob[18];
Lattice lat(N, N);
int M_total, E_total;

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
int heat_capacity() {
	int HC = 0;
	return HC;
}
int magnetic_autocorrelation() {
	int MAC = 0;
	return MAC;
}

void set_probabilities() {
	for (int i = 2; i < 18; i += 2) {
		prob[i] = exp(-beta * i);
	}
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

void sweep_old() {
	int i;
	int nn, sum, delta;
	for (int k = 0; k < N; k++) {
		// Choose a site
		i = std::uniform_int_distribution<int>(0, N * N)(rng);
		// Calculate the sum of the neighbouring spins
		if ((nn = i + 1) >= N) nn -= N;
		sum = lat.at(nn);
		if ((nn = i - 1) < 0) nn += N;
		sum += lat.at(nn);
		if ((nn = i + N) >= N) nn -= N;
		sum += lat.at(nn);
		if ((nn = i - N) < 0) nn += N;
		sum += lat.at(nn);
		// Calculate the change in energy
		delta = sum * lat.at(i);
		//Decide whether to flip spin
		if (delta <= 0) {
			lat.set_at(i, -lat.at(i));
		}
		else if (random_zero_one(rng) < prob[delta]) {
			lat.set_at(i, -lat.at(i));
		}
	}
}

void flip(int x, int y) {
	lat.set(x, y, -lat(x, y));
}



void sweep() {
	for (int i = 0; i < lattice_size; i++) {
		int random_x = std::uniform_int_distribution<int>(0, N - 1)(rng);
		int random_y = std::uniform_int_distribution<int>(0, N - 1)(rng);
		int sum = 0;
		//for (int n : {i - N - 1, i - N, i - N + 1, i - 1, i + 1, i + N - 1, i + N, i + N + 1}) {}
		//above could be faster, but checking periodic boundary cond's becomes messy.
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (!(dx == 0 && dy == 0)) {
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
		}
		int delta_E = 2 * J * lat(random_x, random_y) * sum;
		if (delta_E <= 0) {
			flip(random_x, random_y);
			E_total += delta_E;
			M_total += 2 * lat(random_x, random_y);
		}
		else if (random_zero_one(rng) < prob[delta_E]) {
			flip(random_x, random_y);
			E_total += delta_E;
			M_total += 2 * lat(random_x, random_y);
		}
	}
}

int main() {
	const int mc_steps = 300000;
	const int output_steps = 1000;
	lat.init(1);
	//init_infinite_temperature();
	temperature = 4.0;
	beta = 1. / temperature;
	J = 1;
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	set_probabilities();
	M_total = magnetisation();
	E_total = energy();
	std::ofstream measure_file;
	measure_file.open("measurements_M_E.output");
	measure_file << M_total << ',' << E_total << "\n";
	for (int i = 0; i < mc_steps; ++i) {
		sweep();
		measure_file << M_total << ',' << E_total << "\n";
		if(i%10==0) write_data(i);
		if (i % output_steps == 0) {			
			cout << "m = " << (double)M_total/lattice_size << ", E = " << (double)E_total/lattice_size << ", \n";
		}
		if (i % 10000 == 0) {
			temperature += 0.1;
			beta = 1. / temperature;
			set_probabilities();
			cout << "t = " << temperature << ", " << (char)225 << " = " << beta << "\n";
		}
	}
	measure_file.close();
	return 0;
}