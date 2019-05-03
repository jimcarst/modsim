// (c) Jim Carstens 2019
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
using std::cout; using std::endl;
#pragma warning(disable : 4996)
#include <stdio.h>
class Lattice {
public:
private:
	int _lattice[];
};
const int N = 50;
int lattice[N * N] = { 0 };

void print() {
	for (int i = 0; i < 10; i++) {
		cout << "---";
	}
	cout << "\n";
	for (int i = 0; i < 100; i++) {
		cout << std::showpos << lattice[i] << ", ";
		if ((i + 1) % 10 == 0) {
			cout << "\n";
		}
	}
	for (int i = 0; i < 10; i++) {
		cout << "---";
	}
}

void writeConfiguration(int step) {
	char buffer[128];
	sprintf(buffer, "configurations/conf_step%d.csv", step);
	FILE* fp = fopen(buffer, "w");
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < N - 1; l++) {
			fprintf(fp, "%d,", lattice[k * N + l]);
	
		}
		fprintf(fp, "%d\n", lattice[k * N + N - 1]);
	}
	fclose(fp);
}

int magnetisation() {
	int M = 0;
	for (int i = 0; i < N * N; i++) {
		M += lattice[i];
	}
	return M;
}
int main() {
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<double> random_zero_one(0.0, 1.0);
	const int mc_steps = 1000;

	for (int i = 0; i < N * N; i++) {
		if (i % 2 == 0) {
			lattice[i] = -1;
		}
		else {
			lattice[i] = 1;
		}
	}

	print();
	for (int i = 0; i < mc_steps; ++i) {
		if (random_zero_one(rng) < 0.1) {
			int random_site = std::uniform_int_distribution<int>(0, N * N)(rng);
			lattice[random_site] = -1;
		}
		//print();
		cout << magnetisation() << ", \n";
		writeConfiguration(i);
		//write_data(i);
	}


	return 0;
}