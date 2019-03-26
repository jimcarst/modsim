#include "SudokuSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>

SudokuSolver::SudokuSolver() {
	_dim = 9;
	_sudoku.resize(_dim, vector<int>(_dim, 0));
}


SudokuSolver::~SudokuSolver()
{
}

void SudokuSolver::print() {
	vector<vector<int>> vec = _sudoku;
	cout << setw(28) << setfill('=') << " " << endl<<setfill(' ');
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec[i].size(); j++) {
			if (j % 3 == 0) cout << '|';
			cout << setw(3) << left << vec[j][i];
		}
		if ((i + 1) % 3 == 0) {
			cout << endl << setw(28) << setfill('-') << " " << endl << setfill(' ');
		}
		else {
			cout << endl;
		}
	}
}

void SudokuSolver::read(const char* filename) {
	ifstream ifs(filename);
	if (ifs.fail()) {
		cout << "File doesn't exist" << endl;
		return;
	}
	int value, count(0);
	while (ifs >> value) {
		int x = count % 9;
		int y = count / 9;
		_sudoku[x][y] = value;
		count++;
	}
}

void SudokuSolver::randomChange() {
	for (unsigned int i = 0; i < _sudoku.size(); i++) {
		for (unsigned int j = 0; j < _sudoku[i].size(); j++) {
			_sudoku[i][3] = 9;
		}
	}
}

vector<vector<int>> SudokuSolver::getSu() {
	return this->_sudoku;
}
