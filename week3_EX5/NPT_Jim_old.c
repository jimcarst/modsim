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

double energy() {
	return 0.0;
}
int change_volume(void) {
	double eOld = energy();
	double vOld = box[0] * box[1] * box[2];
	double vNew = vOld + (dsfmt_genrand() * 2 - 1)*deltaV;
	if (vNew < 0) {
		return 0;
	}
	double eNew = energy();
	double arg = -beta * ((eNew - eOld)
		+ P * (vNew - vOld)
		- n_particles * log(vNew / vOld) / beta);
	if (dsfmt_genrand() > exp(arg)) {
		return 0;
	}
	double boxn[NDIM];
	for (int d = 0; d < NDIM; d++) {
		boxn[d] = pow(vNew, 1 / 3.);
	}
	double rtemp[N][NDIM];
	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; d++) {
			rtemp[i][d] = r[i][d] * boxn[d] / box[d];
		}
	}
	/*
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < i; j++) {
			double distance = 0;
			for (int d = 0; d < NDIM; d++) {
				distance += (rtemp[j][d] - rtemp[i][d])*(rtemp[j][d] - rtemp[i][d]);
			}
			if (distance < diameter && i != j) {
				return 0;
			}
		}
	}*/

	for (int i = 0; i < n_particles; i++) {
		for (int d = 0; d < NDIM; d++) {
			r[i][d] = rtemp[i][d];
		}
	}
	for (int d = 0; d < NDIM; d++) {
		box[d] = boxn[d];
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
	for (int i = 0; i < 20; i++) {
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

		set_packing_fraction();

		dsfmt_seed((unsigned int)time(NULL));

	//	printf("#Step \t Volume \t Move-acceptance\t Volume-acceptance \n");

		int move_accepted = 0;
		int vol_accepted = 0;
		int step, n;
		double volumeSum = 0;
		for (step = 0; step < mc_steps; ++step) {
			for (n = 0; n < n_particles; ++n) {
				move_accepted += move_particle();
			}
			vol_accepted += change_volume();

			double thisVolume = box[0] * box[1] * box[2];
			volumeSum += thisVolume;
		/*
			if (step % output_steps == 0) {
				printf("%d \t %lf \t %lf \t %lf \n",
					step, box[0] * box[1] * box[2],
					(double)move_accepted / ((double)n_particles * output_steps),
					(double)vol_accepted / output_steps);
				move_accepted = 0;
				vol_accepted = 0;
				write_data(step);
			}
			*/
		}
		write_data(step);
		double meanVolume = volumeSum / mc_steps;
		printf("{%lg, %lg}, \n", P, meanVolume);
		P += 0.5;
	}
	return 0;
}
