// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>

using std::cout; using std::endl;

int main() {
	int mc_steps = 500000;
	int output_steps = 10000;
	double TStart = 10.;
	double alpha = 0.9;

	SudokuSolver su(9);
	su.read("sudoku.dat");
	cout << "The following sudoku was read in: " << endl;
	su.print();
	su.fillRandom();
	cout << "The sudoku is randomly filled and looks as follows: " << endl;
	su.print();
	int E = su.calculateEnergy();
	su.setTemperature(TStart);
	double temperature = su.getTemperature();
	cout << "E0 = " << E << ", T0 = " << temperature << endl;
	while (E != 0) {
		for (int i = 0; i < mc_steps; i++) {
			su.randomChange();
			E = su.calculateEnergy();
			//su.print();			
			if (i% output_steps == 0 && i != 0) {
				temperature *= alpha;
				su.setTemperature(temperature);
			}
			if (i % 1000 == 0) {
				cout << "{" << temperature << ", " << E << "}, " << endl;
			}			
		}
		su.fillRandom();
		su.setTemperature(TStart);
		temperature = su.getTemperature();
	}
	cout<<"E = " << E <<", T =  "<< temperature<<endl;
	su.print();
	
	return 0;
}