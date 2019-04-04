#include <stdio.h>
#include <time.h>
#include "mt19937.h"

int main() {
	dsfmt_seed(time(NULL));
	printf("random number %lf \n", dsfmt_genrand());
	printf("HOI HOI \n");

	int isPrime(int number) {
		for (int i = 2; i < number; i++) {
			if (number % i == 0) {
				return 0;
			}
		}
		return 1;
	}


	for (int i = 0; i < 10; i++) {
		if (isPrime(i)) {
			printf("%d ", i);
		}
	}

	return 0;
}