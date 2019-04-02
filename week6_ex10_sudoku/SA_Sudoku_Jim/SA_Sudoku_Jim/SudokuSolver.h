#pragma once
#include <vector>
#include <random>
using std::vector;
using std::mt19937;


class SudokuSolver
{
public:
	SudokuSolver(int dim);
	~SudokuSolver();


	void print();
	void read(const char* filename);

	void fillRandom();
	void randomChange();

	int calculateEnergy();
	int rowUniques();
	int colUniques();

	vector<int> blockMaker(int No);
	int blockUniques();
	bool isUnique(vector<int> v, int n);
	vector<vector<int>> getSu();
	//save
	void saveToSu();
	void saveToTemp();
	void save();
	//accept

	double returnRandom();

private:
	int blockFinder(int x, int y);

	vector<vector<int>> _sudoku;
	vector<vector<int>> _sudokuTemp;

	vector<vector<bool>> _fixed;
	int _dim;
	mt19937 _rng;

};

