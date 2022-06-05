#include "../CPP/library/myvector.hpp"
#include "../CPP/library/mymatrix.hpp"

#include <iostream>

int main() {
	matdoub mat(2, 2, 1);
	std::cout << mynorm(mat);
	return 0;
}