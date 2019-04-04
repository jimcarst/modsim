//(c) Jim Carstens 2019
// MC of a NVT ensemble
#include <stdio.h>
#include <math.h>


#define NDIM 3
#define N 10000
#define Nx 3
#define l 1.41421356237309504880168872420969808
#define SITE 4

int n_sites = Nx*Nx*Nx;
int n_particles = Nx*Nx*Nx*SITE;
double diameter = 1.0;
double length = Nx*l;

double r[N][NDIM][SITE];
double box[NDIM];

void write_data() {
	printf("%lg, ", 3.0);
	FILE* fp = fopen("xyz.dat", "w");

	fprintf(fp, "%d\n", n_particles);
	for (int d = 0; d < NDIM; d++) {
		fprintf(fp, "%lg %lg\n", 0.0, box[d]);
	}

	for (int n = 0; n < n_sites; ++n) {
		if (r[10][0][1] == 0) {
			for (int d = 0; d < NDIM; ++d) {
				fprintf(fp, "%lg\t", r[n][d][0]);
				fprintf(fp, "%lg\n", diameter);
			}
		}
		else {
			for (int s = 0; s < SITE; s++) {
				for (int d = 0; d < NDIM; ++d) {
					fprintf(fp, "%lg\t", r[n][d][s]);
				}
				fprintf(fp, "%lg\n", diameter);
			}
		}
	}
	fclose(fp);

}

void makePrim() {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Nx; j++) {
			for (int k = 0; k < Nx; k++) {
				r[i + j * Nx + k * Nx*Nx][0][0] = i * l;
				r[i + j * Nx + k * Nx*Nx][1][0] = j * l;
				r[i + j * Nx + k * Nx*Nx][2][0] = k * l;
			}
		}
	}
}

void makeFcc() {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Nx; j++) {
			for (int k = 0; k < Nx; k++) {
				int loc = i + j * Nx + k * Nx*Nx;
				r[loc][0][0] = i * l;
				r[loc][1][0] = j * l;
				r[loc][2][0] = k * l;

				r[loc][0][1] = i * l + l/2;
				r[loc][1][1] = j * l + l/2;
				r[loc][2][1] = k * l;

				r[loc][0][2] = i * l + l/2;
				r[loc][1][2] = j * l;
				r[loc][2][2] = k * l + l/2;

				r[loc][0][3] = i * l;
				r[loc][1][3] = j * l + l/2;
				r[loc][2][3] = k * l + l/2;
			}
		}
	}
}

int main() {
	printf("%lg, ", 1.0);
	for (int i = 0; i < NDIM; i++) {
		box[i] = length;		
	}
	printf("%lg, ", 8.0);

	/*for (int n = 0; n < Nx; ++n) {
		for (int d = 0; d < NDIM; ++d) {
			r[n][d] = n * l;
		}
	}*/
	printf("%lg, ", 9.0);

	//makePrim();
	makeFcc();

	printf("%lg, ", 10.0);
	write_data();
	printf("%lg, ", 2.0);

	return 0;
}