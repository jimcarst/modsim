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


SudokuSolver::SudokuSolver(int dim = 9) : _dim(dim) {
	_sudoku.resize(_dim, vector<int>(_dim, 0));
	_fixed.resize(_dim, vector<bool>(_dim, false));
	mt19937 _rng(std::chrono::steady_clock::now().time_since_epoch().count());
	_temperature = 10;
	_EOld = calculateEnergy();
}


SudokuSolver::~SudokuSolver() { }

void SudokuSolver::print() {
	int width = 20;
	cout << setw(width) << setfill('=') << " " << endl << setfill(' ');
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


double SudokuSolver::getTemperature() {
	return _temperature;
}

void SudokuSolver::setTemperature(double T) {
	_temperature = T;
}

int SudokuSolver::getEOld() {
	return _EOld;
}



void SudokuSolver::fillRandom() {

	for (int i = 0; i < _dim; i++) {
		for (int j = 0; j < _dim; j++) {

			int nBlock = blockFinder(j, i);
			vector<int> blockVec(_dim);
			for (int i = 0; i < _dim; i++) {
				blockVec[i] = blockMaker(nBlock, i);
			}

			while (!_fixed[j][i] && _sudoku[j][i] == 0) {
				int random = uniform_int_distribution<int>(1, 9)(_rng);
				if (!isUnique(blockVec, random)) {
					_sudoku[j][i] = random;
				}
			}
		}
	}
	_EOld = calculateEnergy();
}

int SudokuSolver::randomChange() {
	int nx, ny;
	do {
		nx = uniform_int_distribution<int>(0, _dim - 1)(_rng);
		ny = uniform_int_distribution<int>(0, _dim - 1)(_rng);
	} while (_fixed[ny][nx]);

	int oldValue = _sudoku[ny][nx];

	int nx2, ny2;
	do {
		nx2 = (nx / 3) * 3 + uniform_int_distribution<int>(0, 2)(_rng);
		ny2 = (ny / 3) * 3 + uniform_int_distribution<int>(0, 2)(_rng);
	} while ((nx2 == nx && ny2 == ny) || _fixed[ny2][nx2]);


	_sudoku[ny][nx] = _sudoku[ny2][nx2];
	_sudoku[ny2][nx2] = oldValue;


	int ENew = calculateEnergy();
	double acceptance = ((double)ENew - (double)_EOld) / _temperature;
	double random = uniform_real_distribution<double>(0.0, 1.0)(_rng);
	if (random < exp(-acceptance)) {//remove energ<olde
		_EOld = ENew;
		return 1;//accepted
	}
	else {
		_sudoku[ny2][nx2] = _sudoku[ny][nx]; //revert change
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
		const vector<int>& colVec = _sudoku[col];
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
	vector<int> rowVec(_dim); // Specify maximum dimension       
	for (int row = 0; row < _dim; row++) {
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
	vector<int> blockVec(_dim);
	for (int nBlock = 0; nBlock < _dim; nBlock++) {
		for (int i = 0; i < _dim; i++) {
			blockVec[i] = blockMaker(nBlock, i);
		}
		for (int i = 1; i <= _dim; i++) {
			if (isUnique(blockVec, i)) {
				count++;
			}
		}
	}
	return count;
}

int SudokuSolver::blockMaker(int No, int i) {
	int block;
	int xmin = 3 * (No % 3);
	int ymin = 3 * (No / 3);
	int col, row;
	col = xmin + (i % 3);
	row = ymin + (i / 3);
	block = _sudoku[col][row];
	return block;
}

int SudokuSolver::blockFinder(int col, int row) {
	int rowFirst = (row / 3) * 3;
	int colFirst = (col / 3) * 3;
	return rowFirst + colFirst / 3;
}

bool SudokuSolver::isUnique(const vector<int> & v, int n) {
	int count = 0;
	for (int i = 0; i < _dim; i++) {
		if (v[i] == n) {
			count++;
		}
	}
	if (count == 1) {
		return true;
	}
	else {
		return false;
	}
}