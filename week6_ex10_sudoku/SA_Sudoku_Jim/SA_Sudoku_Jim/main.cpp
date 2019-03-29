// (c) Jim Carstens 2019

#include "SudokuSolver.h"
#include <iostream>
using namespace std;

vector<int> blockMaker2(int No) {
	vector<int> block(9);
	int xmin = 3 * (No % 3);
	int ymin = 3 * (No / 3);
	SudokuSolver su1(9);
	su1.read("sudokutest.dat");
	su1.print();
	vector<vector<int>> testSu = su1.getSu();
	int col = xmin;
	int row = ymin;
	for (int i = 0; i < 9; i++) {
		col = xmin + (i % 3);
		row = ymin + (i / 3);
		block[i] = testSu[col][row];
	}
	return block;
}



int main() {
	int mc_steps = 1000;

	SudokuSolver su(9);
	//su.print();
	//su.read("sudoku.dat");
	//su.print();
	//su.randomChange();
	//su.print();
	su.read("sudokutest.dat");
	//su.print();
	//cout << su.calculateEnergy() << endl;
	//su.read("sudokuSolved.dat");
	su.print();
	cout << su.calculateEnergy() << endl;
	for (int i = 0; i < mc_steps; i++) {
		su.randomChange();
		cout << su.calculateEnergy() << endl;
	}

	//block1
	/*vector<int> block1(9);
	vector<vector<int>> testSu = su.getSu();
	int row = 0;
	int col = 0;
	for (int i = 0; i < 9; i++) {
		block1[i] = testSu[col][row];
		col++;
		if (col > 2) {
			col = 0;
			row++;
		}
	}*/

	//vector<int> block1 = blockMaker2(2);
	vector<int> block1 = su.blockMaker(6);
	for (auto el : block1) {
		cout << el << ",  ";
	}
	cout << endl;
	//cout <<"energy "<< su.calculateEnergy() << endl;
	su.read("sudokuSolved.dat");
	su.print();
	cout <<"energy "<< su.calculateEnergy() << endl;
	return 0;
}