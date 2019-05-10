// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <stdio.h>
using std::cout; using std::endl;

std::mt19937 rng;
std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
double beta, J;
class Lattice {
public:
private:
	int _lattice[];
};
const int N = 50;
int lattice[N * N] = { 0 };
double prob[5];
void print() {
	for (int i = 0; i < N; i++) {
		cout << "---";
	}
	cout << "\n";
	for (int i = 0; i < N * N; i++) {
		cout << std::showpos << lattice[i] << ", ";
		if ((i + 1) % N == 0) {
			cout << "\n";
		}
	}
	for (int i = 0; i < N; i++) {
		cout << "---";
	}
}

void write_data(int step) {
	char buffer[128];
	snprintf(buffer, 128, "configurations/conf_step%d.csv", step);
	std::ofstream file;
	file.open(buffer);
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < N - 1; l++) {
			file << lattice[k * N + l] << ',';
		}
		file << lattice[k * N + N - 1] << "\n";
	}
	file.close();
}


int magnetisation() {
	int M = 0;
	for (int i = 0; i < N * N; i++) {
		M += lattice[i];
	}
	return M;//per latice
}
int energy() {
	int E = 0;
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

void initialize() {
	for (int i = 2; i < 5; i += 2) prob[i] = exp(-2 * beta * i);
}


void sweep() {
	int i;
	int nn, sum, delta;
	for (int k = 0; k < N; k++) {
		// Choose a site 
		i = std::uniform_int_distribution<int>(0, N * N)(rng);
		// Calculate the sum of the neighbouring spins 
		if ((nn = i + 1) >= N) nn -= N;
		sum = lattice[nn];
		if ((nn = i - 1) < 0) nn += N;
		sum += lattice[nn];
		if ((nn = i + N) >= N) nn -= N;
		sum += lattice[nn];
		if ((nn = i - N) < 0) nn += N;
		sum += lattice[nn];
		// Calculate the change in energy 
		delta = sum * lattice[i];
		//Decide whether to flip spin 
		if (delta <= 0) {
			lattice[i] = -lattice[i];
		}
		else if (random_zero_one(rng) < prob[delta]) {
			lattice[i] = -lattice[i];
		}
	}
}


int main() {
	beta = -2.3;
	double B = 0.;
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	const int mc_steps = 100000;
	const int output_steps = 100;

	for (int i = 0; i < N * N; i++) {
		if (random_zero_one(rng) <= 0.5) {
			lattice[i] = -1;
		}
		else {
			lattice[i] = 1;
		}
	}
	print();
	initialize();
	for (int i = 0; i < mc_steps; ++i) {
		sweep();
		if (i % output_steps == 0) {
			cout << magnetisation() << ", \n";
			write_data(i);
		}
	}

	return 0;
}