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
		if (x >= _width || x < 0 || y >= _heigth || y < 0) { cout << "ERROR_xy, " << x << ", " << y; }
		return _array[y * _width + x];
	}
	int& at(int i) {
		if (i < 0 || i >= _size) { cout << "ERROR_at, "; }
		return _array[i];
	}
	void set(int x, int y, int value) {
		if (x >= _width || x<0 || y>_heigth || y < 0) { cout << "ERROR_set, "; }
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
double beta;
int J;
const int N = 50;
const int lattice_size = N * N;
//int lattice[N * N] = { 0 };
double prob[18];
Lattice lat(N, N);
int M_total;

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
		M += lat.at(i);
	}
	return M;//per latice
}
int energy() {
	int E = 0;
for(int x=0;x)
	}
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
void flip(int x, int y) {
	lat.set(x, y, -lat(x, y));
}

void sweep() {
	for (int i = 0; i < lattice_size; i++) {
		int random_x = std::uniform_int_distribution<int>(0, N - 1)(rng);
		int random_y = std::uniform_int_distribution<int>(0, N - 1)(rng);
		int sum = 0;
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
		}
		else if (random_zero_one(rng) < prob[delta_E]) {
			flip(random_x, random_y);
		//	cout << delta_E<<", ";
		}
	}
}


	int main() {
		//lat.init(-1);
		init_infinite_temperature();
		beta = 1 / 0.5;
		J = 1;
		std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
		const int mc_steps = 100000;
		const int output_steps = 100;
		
		//print();
		initialize();
		M_total = magnetisation();
		for (int i = 0; i < mc_steps; ++i) {
			sweep();
				write_data(i);
			if (i % output_steps == 0) {
				cout<<"M = " << magnetisation() << ", \n";
			}
		}

		return 0;
	}