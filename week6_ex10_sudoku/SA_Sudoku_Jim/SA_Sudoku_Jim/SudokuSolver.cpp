#include "SudokuSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>

using std::vector; 
using std::cout; using std::endl;
using std::setw; using std::setfill; using std::left;
using std::ifstream;
using std::mt19937; using std::uniform_int_distribution;




SudokuSolver::SudokuSolver(int dim = 9) : _dim(dim){
	_sudoku.resize(_dim, vector<int>(_dim, 0));	
}


SudokuSolver::~SudokuSolver()
{
}

void SudokuSolver::print() {
	int width = 20;
	cout << setw(width) << setfill('=') << " " << endl<<setfill(' ');
	for (int i = 0; i < _dim; i++) {
		for (int j = 0; j < _dim; j++) {
			if (j % 3 == 0) cout << '|';
			cout << setw(2) << left << _sudoku[j][i];
		}
		if ((i + 1) % 3 == 0) {
			cout << endl << setw(width) << setfill('-') << " " << endl << setfill(' ');
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
	mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

	
	for (int i = 0; i < _dim; i++) {
		for (int j = 0; j < _dim; j++) {
			int random = uniform_int_distribution<int>(1, 9)(rng);
			_sudoku[j][i] = random;
		}
	}
}

// 243 = 9*9*3
int SudokuSolver::calculateEnergy() {
	int energy = 243 - (rowUniques() + colUniques() + blockUniques());
	return energy;
}


int SudokuSolver::colUniques() {
	int count = 0;
	for (int col = 0; col < _dim; col++) {
		vector<int> colVec = _sudoku[col];
		for (int i = 1; i <= _dim; i++) {
			if (isUnique(colVec, i)) {
				count++;
			}
		}
	}
	return count;
}
int SudokuSolver::rowUniques() {
	int count = 0;
	for (int row = 0; row < _dim; row++) {

		vector<int> rowVec(_dim);
		for (int i = 0; i < _dim; i++) {
			rowVec[i] = _sudoku[i][row];
		}
		for (int i = 1; i <= _dim; i++) {
			if (isUnique(rowVec, i)) {
				count++;
			}
		}
	}
	return count;
}

int SudokuSolver::blockUniques() {
	int count = 0;
	for (int nBlock = 0; nBlock < _dim; nBlock++) {
		vector<int> blockVec = blockMaker(nBlock);
		for (int i = 1; i <= _dim; i++) {
			if (isUnique(blockVec, i)) {
				count++;
			}
		}
	}
	return count;
}

vector<int> SudokuSolver::blockMaker(int No) {
	vector<int> block(_dim);
	int xmin = 3 * (No % 3);
	int ymin = 3 * (No / 3);
	int col, row;
	for (int i = 0; i < _dim; i++) {
		col = xmin + (i % 3);
		row = ymin + (i / 3);
		block[i] = _sudoku[col][row];
	}
	return block;
}

bool SudokuSolver::isUnique(vector<int> v, int n) {
	int count = 0;
	for (size_t i = 0; i < v.size(); i++) {
		if (v[i] == n) {
			count++;
		}
	}
	if (count == 1) {
		return true;
	} else {
		return false;
	}
}

vector<vector<int>> SudokuSolver::getSu() {
	return this->_sudoku;
}
