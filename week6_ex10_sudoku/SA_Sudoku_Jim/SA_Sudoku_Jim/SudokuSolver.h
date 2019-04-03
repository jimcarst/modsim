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
	int randomChange();

	int calculateEnergy();
	int rowUniques();
	int colUniques();

	vector<int> blockMaker(int No);
	int blockUniques();
	bool isUnique(const vector<int>& v, int n);
	vector<vector<int>> getSu();


	double getTemperature();
	void setTemperature(double T);

private:
	int blockFinder(int x, int y);

	vector<vector<int>> _sudoku;
	vector<vector<bool>> _fixed;
	mt19937 _rng;

	double _temperature;
	int _dim;
};

