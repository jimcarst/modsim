// (c) Jim Carstens 2019
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void vecPrint(vector<vector<int>> vec) {
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec[i].size(); j++)
			cout << setw(3) << left << vec[i][j];
		cout << endl;
	}
}

int main() {
	int _dim = 9;

	vector<vector<int>> _sudoku;
	_sudoku.resize(_dim, vector<int>(_dim, 0));


	vecPrint(_sudoku);


	return 0;
}