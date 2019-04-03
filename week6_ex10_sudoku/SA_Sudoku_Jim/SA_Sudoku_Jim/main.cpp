// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>

using std::cout; using std::endl;
//using std::exp;

int main() {
	int mc_steps = 300000;
	int output_steps = 10000;
	double temperature = 5.;
	double alpha = 0.5;

	SudokuSolver su(9);
	su.read("sudoku.dat");
	su.print();
	su.fillRandom();
	su.print();
	int oldEnergy = su.calculateEnergy();
	int energy;
	for (int i = 0; i < mc_steps; i++) {
		su.saveToTemp();
		su.randomChange();
		//su.print();
		energy = su.calculateEnergy();

		double acceptance = ((double)energy - (double)oldEnergy) / temperature;
		if (su.returnRandom() > exp(-acceptance)) {//remove energ<olde
			su.saveToSu();
		}
		else {

		}
		if (i%output_steps == 0) {
			temperature *= alpha;
		cout<<"{" << energy <<", "<< temperature<<"}, "<<endl;
		}
		oldEnergy = energy;
	}
	cout<<"E = " << energy <<", T =  "<< temperature<<endl;
	su.print();
	
	return 0;
}