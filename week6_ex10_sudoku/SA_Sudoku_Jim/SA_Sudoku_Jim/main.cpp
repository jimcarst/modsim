// (c) Jim Carstens 2019
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include "SudokuSolver.h"

using namespace std;

vector<vector<int>> _sudoku;

void vecPrint(vector<vector<int>> vec) {
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec[i].size(); j++) {
			if (j % 3 == 0) cout << '|';
			cout << setw(3) << left << vec[j][i];
		}
		if ((i+1) % 3 == 0) cout <<endl<< setw(28) << setfill('=') << " " << endl;
		cout<<setfill(' ') << endl;
	}
}

void sudokuReader(const char* filename) {
	ifstream ifs(filename);
	//is also a boolean
	if (ifs.fail()) {
		cout << "File doesn't exist" << endl;
		return;
	}
	int xx = 0;
	int yy = 0;
	int count = 0;
	int value;
	while (ifs>>value) {
		//ifs >> value;
		cout << setw(3) << value << xx << " " << yy << endl;
		_sudoku[xx][yy] = value;
		count++;
		xx = count % 9;
		yy = count / 9;
	}
}

int main() {
	int _dim = 9;

	//_sudoku.resize(_dim, vector<int>(_dim, 0));
	//vecPrint(_sudoku);
	//sudokuReader("sudoku.dat");
	//vecPrint(_sudoku);

	SudokuSolver su;
	su.print();

	return 0;
}