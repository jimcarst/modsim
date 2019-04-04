// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>
#include <fstream>

using std::cout; using std::endl;
using std::ofstream;

int main() {
	/*int mc_steps = 1000000;
	int output_steps = 10000;
	double TStart = 3.;
	double alpha = 0.90;*/

	int mc_steps = 100000;
	int output_steps = 3000;
	double TStart = 3.;
	double alpha = 0.8;

	SudokuSolver su(9);
	su.read("sudokuExtreme.dat");
	cout << "The following sudoku was read in: " << endl;
	su.print();
	su.fillRandom();
	cout << "The sudoku is randomly filled and looks as follows: " << endl;
	su.print();
	int E = su.calculateEnergy();
	su.setTemperature(TStart);
	double temperature = su.getTemperature();
	cout << "E0 = " << E << ", T0 = " << temperature << endl;


	ofstream myfile;
	myfile.open("energy.dat");

	while (E != 0) {
		for (int i = 0; i < mc_steps; i++) {
			if (su.randomChange() == 1) {
				E = su.getEOld();
			}
			if (E == 0) break;
			if (i % output_steps == 0 && i != 0) {
				temperature *= alpha;
				su.setTemperature(temperature);
			}
			if (i % output_steps == 0 && i != 0) {
				cout << "{" << temperature << ", " << E << "}, " << endl;
				myfile << temperature << ',' << E << endl;
			}
		}
		if (E == 0) break;
		cout << "{" << temperature << ", " << E << "}, " << endl;
		su.fillRandom();
		E = su.getEOld();
		su.setTemperature(TStart);
		temperature = su.getTemperature();
	}
	myfile << temperature << ',' << E << endl;
	cout << "E = " << E << ", T =  " << temperature << endl;
	su.print();
	myfile.close();
	return 0;
}