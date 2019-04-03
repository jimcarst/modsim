// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>
#include <fstream>

using std::cout; using std::endl;

int main() {
	int mc_steps = 1000000;
	int output_steps = 10000;
	double TStart = 4.;
	double alpha = 0.95;

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


	std::ofstream myfile;
	myfile.open("energy.dat");

	while (E != 0) {
		for (int i = 0; i < mc_steps; i++) {
			if (su.randomChange() == 1) {
				E = su.getEOld();
			}
			//E = su.calculateEnergy();
			if (E == 0) break;
			//su.print();			
			if (i% output_steps == 0 && i != 0) {
				temperature *= alpha;
				su.setTemperature(temperature);
			}
			if (i % 10000 == 0) {
				cout << "{" << temperature << ", " << E << "}, " << endl;
				myfile << temperature << ',' << E << endl;
			}			
		}
		if (E == 0) break;
		su.fillRandom();
		su.setTemperature(TStart);
		temperature = su.getTemperature();
	}
	myfile << temperature << ',' << E << endl;
	cout<<"E = " << E <<", T =  "<< temperature<<endl;
	su.print();
	myfile.close();
	return 0;
}