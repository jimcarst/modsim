#pragma once
#include <iostream>
using std::cout; using std::endl;
class Lattice {
public:
	Lattice(unsigned width, unsigned heigth) : _width(width), _heigth(heigth), _array(0) {
		_size = width * heigth;
		if (width > 0 && heigth > 0) {
			_array = new int[_size];
		}
		init(1);
	}
	~Lattice() {
		delete[] _array;
	}
	void init(int value) {
		for (unsigned i = 0; i < _size; i++) {
			_array[i] = value;
		}
	}
	const int& operator()(int x, int y) const {
		if (x >= _width || x<0 || y>_heigth || y < 0) { return 999; }
		return _array[y * _width + x];
	}
	int& operator()(int x, int y) {
		if (x >= _width || x < 0 || y >= _heigth || y < 0) { cout << "ERROR_xy, " << x << ", " << y; }
		return _array[y * _width + x];
	}
	int& at(int i) {
		if (i < 0 || i >= _size) { cout << "ERROR_at, "; }
		return _array[i];
	}
	void set(int x, int y, int value) {
		if (x >= _width || x<0 || y>_heigth || y < 0) { cout << "ERROR_set, "; }
		_array[y * _width + x] = value;
	}
	void set_at(int i, int value) {
		if (i < 0 || i >= _size) { cout << "ERROR_set_at, "; }
		_array[i] = value;
	}
	// get dims
	unsigned get_width() const { return _width; }
	unsigned get_heigth() const { return _heigth; }

	// private data members
private:
	unsigned _width, _heigth, _size;
	int* _array;
	// to prevent unwanted copying:
	Lattice(const Lattice&);
	Lattice& operator= (const Lattice&);
};