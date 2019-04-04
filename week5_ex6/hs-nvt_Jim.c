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
#define N 520
#define NBINS 400

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 1000;
double packing_fraction = 0.05;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "fcc_108.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];

double histogram[NBINS];
double histogramAverage[NBINS];


void prAr(double arr[], int n) {
	printf("{");
	for (int i = 0; i < n - 1; i++) {
		printf("%lg, ", arr[i]);
	}
	printf("%lg}\n", arr[n - 1]);
}

void initialiseAr(double arr[], int n) {
	for (int i = 0; i < n; i++) {
		arr[i] = 0.;
	}
}

void makeHistogram() {
	double deltaR = 1.0*box[0] / ((double)NBINS);
	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < i; j++) {//i or n_particles
			double r_ij = 0.0;
			if (i == j) continue;
			for (int d = 0; d < NDIM; d++) {
				double r_ij_d = r[i][d] - r[j][d];
				//r_ij_d -= (int)(2.0 * r_ij_d / box[d]) * box[d];
				if (r_ij_d > 0.5*box[d]) {
					r_ij_d -= box[d];
				}
				if (r_ij_d > 0.5*box[d]) {
					double a = 1.0;
				}
				if (r_ij_d < -0.5*box[d]) {
					r_ij_d += box[d];
				}
				assert(r_ij_d <= 0.5 * box[d]);
				assert(r_ij_d >= -0.5 * box[d]);
				r_ij += r_ij_d * r_ij_d;
			}
			r_ij = sqrt(r_ij);
			int nBin = (int)(r_ij / deltaR);
			if (nBin < NBINS) {
				histogram[nBin] += 2.;
			}
			else {
				printf("out of bounds %d", nBin);
			}
		}
	}

	double volume = box[0] * box[1] * box[2];
	double rho = n_particles / volume;
	double preFactor = (4.*M_PI*rho) / 3.;
	for (int i = 0; i < NBINS; i++) {
		double rLower = i * deltaR;
		double rUpper = rLower + deltaR;
		double nIdeal = preFactor * (pow(rUpper, 3.) - pow(rLower, 3.));
		histogram[i] /= (n_particles * nIdeal);
	}

}



/* Read the initial configuration provided */
void read_data(void) {
	FILE* fp = fopen(init_filename, "r");
	if (fp == NULL) {
		printf("Error! File doesn't exist \n");
		return;
	}
	double dmin, dmax, dia;
	fscanf(fp, "%d", &n_particles);
	assert(n_particles < 10000);
	for (int d = 0; d < NDIM; ++d) {
		fscanf(fp, "%lg %lg", &dmin, &dmax);
		box[d] = fabs(dmax - dmin);
	}
	for (int n = 0; n < n_particles; ++n) {
		for (int d = 0; d < NDIM; ++d) fscanf(fp, "%lg ", &r[n][d]);
		fscanf(fp, "%lg", &dia);
		assert(dia == diameter);
	}
	fclose(fp);
}

/* Move a particle randomly */
int move_particle(void) {
	//Choose a random particle
	int rpid = (unsigned int)(dsfmt_genrand()*n_particles);
	if (rpid >= n_particles || rpid < 0) {//here???
		return 0;
	}

	double new_pos[NDIM];
	int d;
	for (d = 0; d < NDIM; ++d) {
		//Displace by a random value between -delta and delta
		new_pos[d] = r[rpid][d] + delta * (2.0 * dsfmt_genrand() - 1.0);
		//Apply periodic boundaries
		//new_pos[d] -= floor(new_pos[d] / box[d]) * box[d];
		if (new_pos[d] > box[d]) {
			new_pos[d] -= box[d];
		}
		if (new_pos[d] < 0.0) {
			new_pos[d] += box[d];
		}
		assert(new_pos[d] < box[d]);
		assert(new_pos[d] >= 0.0);
	}

	int n;
	//Check for overlaps
	for (n = 0; n < n_particles; ++n) {
		if (n == rpid) continue;
		double dist = 0.0;
		for (d = 0; d < NDIM; ++d) {
			double min_d = new_pos[d] - r[n][d];
			// Find the distance with the Nearest Image Convension
			//min_d -= (int)(2.0 * min_d / box[d]) * box[d];
			if (min_d > 0.5*box[d]) {
				min_d -= box[d];
			}
			if (min_d < -0.5*box[d]) {
				min_d += box[d];
			}
			assert(min_d <= 0.5 * box[d]);
			assert(min_d >= -0.5 * box[d]);
			dist += min_d * min_d;
		}
		if (dist <= diameter * diameter) {
			//reject the move
			return 0;
		}
	}

	//Accept the move if reached here
	for (d = 0; d < NDIM; ++d) r[rpid][d] = new_pos[d];

	return 1;
}

/* Write the configuration files */
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

/* Scales the box volume and coordinates read from 'init_fileaname'
 * according to the value of 'packing fraction'
 * */
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
	for (int q = 0; q < 10; q++) {
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

		//Read the input configuration
		read_data();

		if (n_particles == 0) {
			printf("Error: Number of particles, n_particles = 0.\n");
			return 0;
		}

		initialiseAr(histogramAverage, NBINS);

		//Set the packing fraction according to the set variable
		set_packing_fraction();
		//printf("Box: %lf %lf %lf\n", box[0], box[1], box[2]);

		dsfmt_seed((unsigned int)time(NULL));

		int accepted = 0;
		int step, n;
		//Perform MC moves
		for (step = 0; step < mc_steps; ++step) {
			for (n = 0; n < n_particles; ++n) {
				accepted += move_particle();
			}

			if (step % output_steps == 0) {
				//printf("Step %d. Move acceptance: %f.\n", step, (double)accepted / ((double)n_particles * output_steps));
				accepted = 0;
				write_data(step);
			}
			initialiseAr(histogram, NBINS);
			makeHistogram();
			for (int i = 0; i < NBINS; i++) {
				histogramAverage[i] += histogram[i] / ((double)mc_steps);
			}
		}

		/*char buffer[128];
		sprintf(buffer, "measurement_pack%03d.output", (int)(packing_fraction*100.));
		FILE* fp = fopen(buffer, "w");
		for (int pi = 0; pi < NBINS; pi++) {
			fprintf(fp, "%lg\n", histogramAverage[pi]);
		}
		fclose(fp);*/
		printf("pack %lg length %lg\n", packing_fraction, box[0]);
		//prAr(histogramAverage, NBINS);

		packing_fraction += 0.05;

	}
	getchar();
	return 0;
}
