#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937.h"
#include <windows.h>
#include <math.h>

//Markov chain method lattice

#define M 10

void prAr(double arr[], int n) {
	int len = n * n;
	printf("{ {");
	for (int i = 0; i < len - 1; i++) {
		if ((i + 1) % n != 0) {
			printf("%lg, ", arr[i]);
		} else {
			printf("%lg},\n{", arr[i]);
		}
	}
	printf("%lg} }\n", arr[len - 1]);
}
double mean(double arr[], int len) {
	double sum = 0;
	for (int i = 0; i < len; i++) {
		sum += arr[i];
	}

	return sum / len;
}

double variance(double arr[], int len) {
	double sd = 0;
	double mu = mean(arr,len);
	for (int i = 0; i < len; i++) {
		sd += (arr[i] - mu)*(arr[i] - mu);
	}
	sd /= (len - 1);
	return sd;
}


int nb(int len, int site, int k) {
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

double amountnb(int len, int site) {
	if (site == 0
		|| site == len - 1
		|| site == len * len - 1
		|| site == len * (len - 1)) {
		return 2.;
	}
	else if (((site + 1) % len == 0)
		|| !(site + 1 <= len * (len - 1))
		|| site % len == 0
		|| !(site + 1 >= len)) {
		return 3.;
	}
	else {
		return 4.;
	}
}
int main() {
	int dim = M * M;
	double grid[M*M];
	for (int i = 0; i < dim; i++) {
		grid[i] = 0;
	}
	int position = dim - 1;
	grid[position] = 1;
	//prAr(grid, M);	
	dsfmt_seed(time(NULL));

	double no;
	int newPosition;
	int ntrial = 10;
	for (int i = 0; i < ntrial; i++) {
		no = floor(dsfmt_genrand() * 4 + 1);
		int neighbour = nb(M, position, no);
		if (neighbour == -1) {
			newPosition = position;
		}
		else {
			newPosition = neighbour;
		}
		grid[newPosition]++;
		position = newPosition;
		//double stdev = sqrt(variance(grid, dim)) / mean(grid, dim);
		//printf("%lg, ", stdev);
	}
		//prAr(grid, M);
		
		double tM[M*M*M*M];
		for (int i = 0; i < M*M*M*M; i++) {
			tM[i] = 0;
		}
		for (int i = 0; i < M*M*M*M; i++) {
			int from = i % (M*M);
			int to = i / (M*M);
			for (int k = 1; k < 5; k++) {
				if (nb(M, from, k) == to) {
					tM[i] += 0.25;
				}
			}
			if (from == to) {
				tM[i] = 1 - amountnb(M, from) / 4;
			}
		}/*
		for (int i = 0; i < M*M*M*M; i++) {
			int from = i % (M*M);
			int to = i / (M*M);
			tM[i] = amountnb(M, from);
		}*/
		prAr(tM, M*M);
		
		
	return 0;
}