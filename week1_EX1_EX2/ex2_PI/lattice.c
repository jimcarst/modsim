#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937.h"
#include <windows.h>
#include <math.h>

//Markov chain method lattice

#define M 3

void prAr(double arr[], int n) {
	int len = n * n;
	for (int i = 0; i < len - 1; i++) {
		printf("%lg, ", arr[i]);
		if ((i + 1) % n == 0) {
		printf("\n");
		}
	}
	printf("%lg \n", arr[len - 1]);
}
void nbprint(int len, int site) {
	printf("{%d: ", site);
	//1
	if ((site + 1) % len == 0 && site != 0) {
		printf("%d, ", -1);
	} else {
		printf("%d, ", site + 1);
	}
	//2
	if (site + 1 > len*(len - 1)) {
		printf("%d, ", -1);
	} else {
		printf("%d, ", site + len);
	}
	//3
	if (site % len == 0 || site == 0) {
		printf("%d, ", -1);
	} else {
		printf("%d, ", site - 1);
	}
	//4
	if (site + 1 < len) {
		printf("%d", -1);
	} else {
		printf("%d", site- len);
	}
	printf("},\n");
}

int nb2(int len, int site, int k) {
	if (k == 1 && !((site + 1) % len == 0 && site != 0)) {
		return site + 1;
	}
	else if (k == 2 && (site + 1 <= len * (len - 1))) {
		return site + len;
	}
	else if (k == 3 && !(site % len == 0 || site == 0)) {
		return site - 1;
	}
	else if (k == 4 && site + 1 >= len) {
		return site - len;
	}
	else return -1;
}

int amountnb(int len, int site) {
	if (site == 0
		|| site == len - 1 
		|| site == len * len - 1 
		|| site == len * (len - 1) - 1) {
		return 2;
	}
	else if (((site + 1) % len == 0)
		|| !(site + 1 <= len * (len - 1))
		|| site % len == 0
		|| !(site + 1 >= len)) {
		return 3;
	}
	else {
		return 4;
	}
}
int main() {
	int dim = M * M;
	double grid[M*M];
	for (int i = 0; i < dim; i++) {
		grid[i] = 0;
	}
	grid[dim - 1] = 1;
	prAr(grid, M);
	for (int i = 0; i < dim; i++) {
		nbprint(M, i);
	}
	for (int j = 0; j < dim; j++) {
		for (int i = 1; i < 5; i++) {
			printf(",  %d", nb2(M, j, i));
		}
	printf("\n");
	}

	prAr(grid, M);
	
	for (int i = 0; i < dim; i++) {
		double nNeighbours = amountnb(M, i);
		for (int k = 1; k < 5; k++) {
			int neighbour = nb2(M, i, k);
			double valueold = grid[i];
			if (neighbour != -1) {
				grid[neighbour] += valueold / 4;
				
			}
		}
		grid[i] = grid[i] - grid[i] / nNeighbours;
	}
	prAr(grid, M);

	////  	
	double tM[M*M*M*M];
	for (int i = 0; i < M*M*M*M; i++) {
		int x = i % (M*M);
		int y = i / (M*M);
		if (abs(x - y) > 3) {
			tM[i] = 0;
		}
		else {
			tM[i] = 1;
		}
	}
	//prAr(tM, M*M);

	double no;

	dsfmt_seed(time(NULL));
	no = floor(dsfmt_genrand() * 4 + 1);
	printf("%lg", no);

	return 0;
}