//#include <time.h>
//#include <assert.h>
//#include "mt19937.h"
#include <stdio.h>
#include <math.h>
#include "sinfft.h"
#pragma warning(disable : 4996)


#define M_PI 3.14159265358979323846
#define N 2048
#define length 10.
#define ETA 0.4

double g_r_Old[N];
double c_r[N];
double c_q[N];
double g_q[N];
double g_r_New[N];
double S_q[N];

double sigma = 1.;
double deltaR = length / N;
double deltaQ = M_PI / N / (length / N);
double rho = 6. / M_PI * ETA; //set a value!
double TOL = 0.0001;
int maxIterations = 100000;
double alpha = 0.5;

void prAr(double arr[], int n) {
	printf("{");
	for (int i = 0; i < n - 1; i++) {
		printf("%lg, ", arr[i]);
	}
	printf("%lg}\n", arr[n - 1]);
}
void initAr(double arr[], int len, double val) {
	for (int i = 0; i < len; i++) {
		arr[i] = val;
	}
}

void copyAr(double from[], double to[], int len) {
	for (int i = 0; i < len; i++) {
		to[i] = from[i];
	}
}

void initAll() {
	initAr(g_r_Old, N, 0.);
	initAr(c_r, N, 0.);
	initAr(c_q, N, 0.);
	initAr(g_q, N, 0.);
	initAr(g_r_New, N, 0.);
}

void compute_c() {
	for (int i = 0; i < N; i++) {
		if (i*deltaR < sigma) {
			c_r[i] = -g_r_Old[i] - 1.;	//??maybe g_r_Old
			//c_r[i] = -g_r_Old[i] - deltaR*i;
		}
		else {
			c_r[i] = 0.;
		}
	}
}

void compute_g_q() {
	for (int i = 0; i < N; i++) {
		//g_q[i] = (rho * c_q[i] * c_q[i]) / (1. - rho * c_q[i]);
		//g_q[i] = (rho * c_q[i] * c_q[i]) / ((double)i*deltaR -rho * c_q[i]);
		g_q[i] = deltaQ * deltaQ / (2.0*M_PI*M_PI*deltaR) / i * rho*c_q[i] * c_q[i] / (1.0 - rho * c_q[i] / i);
	}
}

double compute_difference() {
	double diff = 0.;
	for (int i = 0; i < N; i++) {
		diff += fabs(g_r_New[i] - g_r_Old[i]);
	}
	return diff;
}

void OZstep(int print) {

	if (print == 1) {
		printf("g_r_Old\n");
		prAr(g_r_Old, N);
	}
	compute_c();
	if (print == 1) {
		printf("c_r\n");
		prAr(c_r, N);
	}
	copyAr(c_r, c_q, N);
	//edit c_q
	for (int i = 0; i < N; i++) {
		c_q[i] = i * 4 * M_PI *deltaR*deltaR / deltaQ * c_q[i];
	}
	c_q[0] = 0;
	sinft(c_q, N);
	if (print == 1) {
		printf("FFT(c_r)\n");
		prAr(c_q, N);
	}

	compute_g_q();
	if (print == 1) {
		printf("g_q\n");
		prAr(g_q, N);
	}

	copyAr(g_q, g_r_New, N);
	g_q[0] = 0;
	sinft(g_r_New, N);
	for (int i = 0; i < N; i++) {
		g_r_New[i] = g_r_New[i] / i;
	}
	if (print == 1) {
		printf("FFT(g_q)\n");
		prAr(g_r_New, N);
	}
}
void write_data_g(int step) {
	char buffer[128];
	sprintf(buffer, "meas_g_pack%03d.output", step);
	FILE* fp = fopen(buffer, "w");
	for (int i = 0; i < N; i++) {
		fprintf(fp, "%lg\n", g_r_Old[i]);
	}
	fclose(fp);
}
void write_data_c(int step) {
	char buffer[128];
	sprintf(buffer, "meas_c_pack%03d.output", step);
	FILE* fp = fopen(buffer, "w");
	for (int i = 0; i < N; i++) {
		fprintf(fp, "%lg\n", c_r[i]);
	}
	fclose(fp);
}
void write_data_S(int step) {
	for (int i = 0; i < N; i++) {
		S_q[i] = 1.0 / (1.0 - rho * c_q[i] / i);
	}

	char buffer[128];
	sprintf(buffer, "meas_S_pack%03d.output", step);
	FILE* fp = fopen(buffer, "w");
	for (int i = 0; i < N; i++) {
		fprintf(fp, "%lg\n", S_q[i]);
	}
	fclose(fp);
}

int main(int argc, char* argv[]) {

	initAll();
	int done = 0;
	int count = 0;
	while (done == 0 && count < maxIterations) {
		OZstep(0);
		if (compute_difference() <= TOL) {
			printf("done");
			done = 1;
			for (int i = 0; i < N; i++) {
				g_r_Old[i] = g_r_New[i];
			}
			prAr(g_r_Old, N);
		}
		else {
			//printf("not done!");
			done = 0;
			//copyAr(g_r_New, g_r_Old, N);
			for (int i = 0; i < N; i++) {
				g_r_Old[i] = alpha * g_r_New[i] + (1 - alpha)*g_r_Old[i];
			}
		}
		count++;
	}
	printf("rho: %lf, gp: %d, dr: %lf\n", rho, N, deltaR);
	prAr(g_r_Old, N);
	prAr(c_r, N);
	write_data_g((int)(ETA * 100));	
	write_data_c((int)(ETA * 100));
	write_data_S((int)(ETA * 100));
    return 0;
}
