#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#pragma warning(disable : 4996)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 10000;
const int output_steps = 100;
const double packing_fraction = 0.2;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];


/* Functions */
void prAr(double arr[], int n) {
	printf("{");
	for (int i = 0; i < n - 1; i++) {
		printf("%lg, ", arr[i]);
	}
	printf("%lg}\n", arr[n - 1]);
}
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

void read_data(void) {
	FILE *fp = fopen(init_filename, "r");
	if (fp == NULL) {
		printf("Error! File doesn't exist \n");
		return;
	}
	int npart;
	fscanf(fp, "%d", &npart);
	n_particles = npart;
	for (int d = 0; d < NDIM; d++) {
		double dim;
		fscanf(fp, "%*lg %lg", &dim);
		box[d] = dim;
	}
	for (int i = 0; i < n_particles; i++) {
		double x, y, z;
		fscanf(fp, "%lg %lg %lg %*lg", &x, &y, &z);
		r[i][0] = x;
		r[i][1] = y;
		r[i][2] = z;
	}
}


int move_particle(void) {
	int pos = (int)(dsfmt_genrand()*n_particles);
	double dxi[NDIM];
	double postemp[NDIM];

	for (int d = 0; d < NDIM; d++) {
		dxi[d] = (dsfmt_genrand() * 2 - 1) * delta;
		postemp[d] = r[pos][d] + dxi[d];
		//by using this order, the periodic boundary conditions
		//automatically apply for the overlap check
		if (postemp[d] < 0) {
			postemp[d] += box[d];
		}
		if (r[pos][d] + dxi[d] >= box[d]) {
			postemp[d] -= box[d];
		}
	}
	for (int i = 0; i < n_particles; i++) {
		double distance = 0;
		for (int d = 0; d < NDIM; d++) {
			distance += (postemp[d] - r[i][d])*(postemp[d] - r[i][d]);
		}
		if (distance < diameter && i != pos) {
			return 0;
		}
	}

	for (int d = 0; d < NDIM; d++) {
		r[pos][d] = postemp[d];
	}
	return 1;
}

void write_data(int step) {
	char buffer[128];
	sprintf(buffer, "coords_step%07d.output", step);
	FILE* fp = fopen(buffer, "w");
	int d, n;
	fprintf(fp, "%d\n", n_particles);
	for (d = 0; d < NDIM; ++d) {
		fprintf(fp, "%lf %lf\n", 0.0, box[d]);
	}
	for (n = 0; n < n_particles; ++n) {
		for (d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
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

	assert(packing_fraction > 0.0 && packing_fraction < 1.0);
	assert(diameter > 0.0);
	assert(delta > 0.0);

	radius = 0.5 * diameter;

	if (NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
	else if (NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
	else {
		printf("Number of dimensions NDIM = %d, not supported.", NDIM);
		return 0;
	}
	//prAr(box, NDIM);
	read_data();
	//prAr(box, NDIM);
	//prAr2D(r, n_particles*NDIM , NDIM);

	if (n_particles == 0) {
		printf("Error: Number of particles, n_particles = 0.\n");
		return 0;
	}

	set_packing_fraction();

	dsfmt_seed(time(NULL));

	int accepted = 0;
	int step, n;
	for (step = 0; step < mc_steps; ++step) {
		for (n = 0; n < n_particles; ++n) {
			accepted += move_particle();
		}

		if (step % output_steps == 0) {
			printf("Step %d. Move acceptance: %lf.\n", step, (double)accepted / (n_particles * output_steps));
			accepted = 0;
			write_data(step);
		}
	}

	return 0;
}
