// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <stdio.h>
using std::cout; using std::endl;

class Lattice {
public:
	Lattice(unsigned width, unsigned heigth) : _width(width), _heigth(heigth), _array(0) {
		_size = width * heigth;
		if (width > 0 && heigth > 0) {
			_array = new int[_size];
		}
		init(1);
	}
	~Lattice() {
		delete[] _array;
	}
	void init(int value) {
		for (unsigned i = 0; i < _size; i++) {
			_array[i] = value;
		}
	}
	const int& operator()(int x, int y) const {
		if (x >= _width || x<0 || y>_heigth || y < 0) { return 999; }
		return _array[y * _width + x];
	}
	int& operator()(int x, int y) {
		if (x >= _width || x < 0 || y >= _heigth || y < 0) { cout << "ERROR"; }
		return _array[y * _width + x];
	}
	int& at(int i) {
		if (i < 0 || i >= _size) { cout << "ERROR"; }
		return _array[i];
	}
	void set(int x, int y, int value) {
		if (x >= _width || x<0 || y>_heigth || y < 0) { cout << "ERROR"; }
		_array[y * _width + x] = value;
	}
	// get dims
	unsigned get_width() const { return _width; }
	unsigned get_heigth() const { return _heigth; }

	// private data members
private:
	unsigned _width, _heigth, _size;
	int* _array;
	// to prevent unwanted copying:
	Lattice(const Lattice&);
	Lattice& operator= (const Lattice&);
};

std::mt19937 rng;
std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
double beta, J;
const int N = 50;
const int lattice_size = N * N;
//int lattice[N * N] = { 0 };
double prob[8];
Lattice lat(N, N);

void print() {
	for (int i = 0; i < N; i++) {
		cout << "---";
	}
	cout << "\n";
	for (int i = 0; i < N * N; i++) {
		cout << std::showpos << lat.at(i) << ", ";
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
			file << lat(k, l) << ',';
		}
		file << lat(k, N - 1) << "\n";
	}
	file.close();
}


int magnetisation() {
	int M = 0;
	for (int i = 0; i < N * N; i++) {
		//	M += lattice[i];
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
	for (int i = 2; i < 5; i += 2) {
		prob[i] = exp(-2 * beta * i);
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
/*
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
*/
void sweep2() {
	for (int i = 0; i < lattice_size; i++) {
		int random_x = std::uniform_int_distribution<int>(0, N)(rng);
		int random_y = std::uniform_int_distribution<int>(0, N)(rng);
		double delta_E = 0.;
		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				if (!(dx == 0 && dy == 0)) {
					delta_E += 2 * J * lat(random_x, random_y) * lat(random_x + dx, random_y + dy);
					//periodic boundary conditions
				}
			}
		}
	}
}


	int main() {
		lat.init(1);
		beta = 1 / 2.;
		std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
		const int mc_steps = 100000;
		const int output_steps = 100;


		print();
		initialize();
		for (int i = 0; i < mc_steps; ++i) {
			sweep2();
			if (i % output_steps == 0) {
				cout << magnetisation() << ", \n";
				write_data(i);
			}
		}

		return 0;
	}