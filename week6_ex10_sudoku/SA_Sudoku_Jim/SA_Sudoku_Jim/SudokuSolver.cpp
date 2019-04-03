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
using std::mt19937; using std::uniform_int_distribution; using std::uniform_real_distribution;




SudokuSolver::SudokuSolver(int dim = 9) : _dim(dim){
	_sudoku.resize(_dim, vector<int>(_dim, 0));	
	_fixed.resize(_dim, vector<bool>(_dim, false));
	mt19937 _rng(std::chrono::steady_clock::now().time_since_epoch().count());
	_temperature = 10;
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

	for (int i = 0; i < _dim; i++) {
		for (int j = 0; j < _dim; j++) {
			if (_sudoku[j][i] != 0) {
				_fixed[j][i] = true;
			}
		}
	}
}

double SudokuSolver::returnRandom() {
	double random = uniform_real_distribution<double>(0.0, 1.0)(_rng);
	return random;
}

double SudokuSolver::getTemperature() {
	return _temperature;
}

void SudokuSolver::setTemperature(double T) {
	_temperature = T;
}

int SudokuSolver::blockFinder(int col, int row) {
	int rowFirst = (row / 3) * 3;
	int colFirst = (col / 3) * 3;
	return rowFirst  + colFirst / 3;
}

void SudokuSolver::fillRandom() {

	for (int i = 0; i < _dim; i++) {
		for (int j = 0; j < _dim; j++) {
			
			int nBlock = blockFinder(j, i);
			vector<int> blockVec = blockMaker(nBlock);

			while (!_fixed[j][i] && _sudoku[j][i] == 0) {
				int random = uniform_int_distribution<int>(1, 9)(_rng);
				if (!isUnique(blockVec, random)) {
					_sudoku[j][i] = random;
				}
			}
		}
	}

}

int SudokuSolver::randomChange() {
	int EOld = calculateEnergy();
	int nx = uniform_int_distribution<int>(0, _dim - 1)(_rng);
	int ny = uniform_int_distribution<int>(0, _dim - 1)(_rng);
	int oldValue = _sudoku[ny][nx];
	if (!_fixed[ny][nx]) {
		_sudoku[ny][nx] = uniform_int_distribution<int>(1, 9)(_rng);
	}

	int ENew = calculateEnergy();
	double acceptance = ((double)ENew - (double)EOld) / _temperature;
	double random = uniform_real_distribution<double>(0.0, 1.0)(_rng);
	if (random < exp(-acceptance)) {//remove energ<olde
		return 1;//accepted
	}
	else {
		_sudoku[ny][nx] = oldValue; //revert change
		return 0; //not accepted
	}

}

// 243 = 9*9*3
int SudokuSolver::calculateEnergy() {
	int energy = 243 - (rowUniques() + colUniques() + blockUniques());//count number as faults
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

void SudokuSolver::saveToSu() {
	_sudoku = _sudokuTemp;
}

void SudokuSolver::saveToTemp() {
	_sudokuTemp = _sudoku;
}

void SudokuSolver::save() {
	_sudoku = _sudokuTemp;
}
