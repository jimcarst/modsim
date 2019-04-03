// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>

using std::cout; using std::endl;

int main() {
	int mc_steps = 300000;
	int output_steps = 10000;
	double temperature = 10.;
	double alpha = 0.8;

	SudokuSolver su(9);
	su.read("sudoku.dat");
	cout << "The following sudoku was read in: " << endl;
	su.print();
	su.fillRandom();
	cout << "The sudoku is randomly filled and looks as follows: " << endl;
	su.print();
	int E = su.calculateEnergy();
	su.setTemperature(temperature);
	cout << "E0 = " << E << ", T0 = " << temperature << endl;
	while (E != 0) {
		for (int i = 0; i < mc_steps; i++) {
			su.randomChange();
			E = su.calculateEnergy();
			//su.print();			
			if (i % output_steps == 0) {
				temperature *= alpha;
				su.setTemperature(temperature);
			}
			if (i % 100 == 0) {
				cout << "{" << E << ", " << temperature << "}, " << endl;
			}
		}
	}
	cout<<"E = " << E <<", T =  "<< temperature<<endl;
	su.print();
	
	return 0;
}