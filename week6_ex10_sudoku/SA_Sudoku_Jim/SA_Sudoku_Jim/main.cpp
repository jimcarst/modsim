// (c) Jim Carstens 2019

#include "SudokuSolver.h"

int main() {

	SudokuSolver su;
	su.print();
	su.read("sudoku.dat");
	su.print();
	su.randomChange();
	su.print();
	return 0;
}