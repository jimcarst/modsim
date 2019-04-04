#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"
#include <assert.h>
#pragma warning(disable : 4996)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;
const double packing_fraction = 0.6;
const double diameter = 1.0;
const double delta = 0.1;
/* Volume change -deltaV, delta V */
const double deltaV = 2.0;
/* Reduced pressure \beta P */
const double betaP = 3.0;
const double beta = 1.0;
double P = 0.1;
const char* init_filename = "fluid2.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];
double overlap[10][2];


void prAr2D(double arr[], int n, int dim) {
	printf("{");
	for (int i = 0; i < n - 1; i++) {
		printf("%lg, ", arr[i]);
		if ((i + 1) % dim == 0) {
			printf("\n");
		}
	}
	printf("%lg}\n", arr[n - 1]);
}

int checkOverlap() {
	int count = 0;
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < i; j++) {
			double distance = 0;
			for (int d = 0; d < NDIM; d++) {
				distance += (r[j][d] - r[i][d])*(r[j][d] - r[i][d]);
			}
			if (distance < diameter && i != j) {
				overlap[count][0] = (double)i;
				overlap[count][1] = (double)j;
				count++;
				return 0;
			}
		}
	}
	return 1;
}

void read_data(void) {
	FILE *fp = fopen(init_filename, "r");
	if (fp == NULL) {
		printf("Error! File doesn't exist \n");
		return;
	}
	int npart;
	assert(fscanf(fp, "%d", &npart) > 0);
	n_particles = npart;
	for (int d = 0; d < NDIM; d++) {
		double dim;
		assert(fscanf(fp, "%*lg %lg", &dim) > 0);
		box[d] = dim;
	}
	for (int i = 0; i < n_particles; i++) {
		double x, y, z;
		assert(fscanf(fp, "%lg %lg %lg %*lg", &x, &y, &z) > 0);
		r[i][0] = x;
		r[i][1] = y;
		r[i][2] = z;
	}
}

void write_data(int step) {
	char buffer[128];
	sprintf(buffer, "coords_step%07d.dat", step);
	FILE* fp = fopen(buffer, "w");
	int d, n;
	fprintf(fp, "%d\n", n_particles);
	for (d = 0; d < NDIM; ++d) {
		fprintf(fp, "%lf %lf\n", 0.0, box[d]);
	}
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
		fprintf(fp, "%lf\n", diameter);
	}
	fclose(fp);
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

int main(int argc, char* argv[]) {
	
		radius = 0.5 * diameter;
		if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
		else if (NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
		else {
			printf("Number of dimensions NDIM = %d, not supported.", NDIM);
			return 0;
		}
		read_data();

		if (n_particles == 0) {
			printf("Error: Number of particles, n_particles = 0.\n");
			return 0;
		}

		for (int i = 0; i < 10; i++) {
			for (int d = 0; d < 2; d++) {
				overlap[i][d] = 0;
			}
		}

		//set_packing_fraction();
		printf("%d\n", checkOverlap());
		prAr2D(*r, n_particles*NDIM, NDIM);
		prAr2D(*overlap, 10*2, 2);
		printf("1, x: %lg, y: %lg, z: %lg\n", r[1][0], r[1][1], r[1][2]);
		printf("2, x: %lg, y: %lg, z: %lg\n", r[35][0], r[35][1], r[35][2]);
	return 0;
}
