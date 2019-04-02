// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>
using namespace std;


int main() {
	int mc_steps = 10000;
	int output_steps = 100;
	double temperature = 5;
	double alpha = 0.99;

	SudokuSolver su(9);
	su.read("sudoku.dat");
	su.print();
	su.fillRandom();
	su.print();
	int oldEnergy = 300;
	int energy;
	for (int i = 0; i < mc_steps; i++) {
		su.saveToTemp();
		su.randomChange();
		//su.print();
		energy = su.calculateEnergy();

		double acceptance = ((double)oldEnergy - (double)energy) / temperature;
		if (energy > oldEnergy || su.returnRandom() > exp(acceptance))	{
			su.saveToSu();
		}
		else {

		}
		if (i%output_steps == 0) {
			temperature *= alpha;
		}
		oldEnergy = energy;
		cout<<"{" << energy <<", "<< temperature<<"}, "<<endl;
	}
	cout<<"E = " << energy <<", T =  "<< temperature<<endl;
	su.print();
	
	return 0;
}