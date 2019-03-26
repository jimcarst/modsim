#pragma once
#include <vector>
using namespace std;

class SudokuSolver
{
public:
	SudokuSolver();
	~SudokuSolver();


	void print();
	void read(const char* filename);
	void randomChange();

	vector<vector<int>> getSu();

private:
	vector<vector<int>> _sudoku;
	int _dim;
};

