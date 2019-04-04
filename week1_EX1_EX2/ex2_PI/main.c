#include <stdio.h>
#include <time.h>
#include "mt19937.h"
#include <windows.h>

//direct sampling methods

int main() {
	for (int loop = 0; loop < 20; loop++) {
		dsfmt_seed(time(NULL));
		int nhit = 0;
		int ntrial = 100;
		for (int i = 1; i < ntrial; i++) {
			double x, y;
			x = dsfmt_genrand() * 2 - 1;
			y = dsfmt_genrand() * 2 - 1;
			if (x * x + y * y < 1) {
				nhit++;
			}
		}
		double result = 4. * nhit / ntrial;
		printf("%lf,", result);
		Sleep(1000);
	}
	return 0;
}