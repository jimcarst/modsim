#pragma once
#include <vector>
using std::vector;

class SudokuSolver
{
public:
	SudokuSolver(int dim);
	~SudokuSolver();


	void print();
	void read(const char* filename);
	void randomChange();

	int calculateEnergy();
	int rowUniques();
	int colUniques();

	vector<int> blockMaker(int No);
	int blockUniques();
	bool isUnique(vector<int> v, int n);
	vector<vector<int>> getSu();
	//save

private:
	vector<vector<int>> _sudoku;
	int _dim;
};

