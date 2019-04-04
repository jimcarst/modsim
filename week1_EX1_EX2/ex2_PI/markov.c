#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937.h"
#include <windows.h>

//Markov chain method

int main() {
	dsfmt_seed(time(NULL));
	int ntrial = 10000;
	double dmax = 0.6;
	for (int loop = 0; loop < 50; loop++) {
		int nhit = 0;
		int accept = 0;
		int reject = 0;
		double x = 0;
		double y = 0;
		for (int i = 1; i < ntrial; i++) {
			double dx, dy;
			dx = (dsfmt_genrand() * 2 - 1)*dmax;
			dy = (dsfmt_genrand() * 2 - 1)*dmax;
			if (abs(x + dx) < 1 && abs(y + dy) < 1) {
				x += dx;
				y += dy;
			}
			if (x * x + y * y < 1) {
				nhit++;
				accept++;
			}
			else
				reject++;
		}
		double result = 4. * nhit / ntrial;
		printf("%lg, ", result);
	//	Sleep(1000);
	}
	return 0;
}