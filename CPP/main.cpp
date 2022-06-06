#include "../CPP/library/myvector.hpp"
#include "../CPP/library/mymatrix.hpp"
#include "../CPP/library/linalg/gaussj.hpp"

int main() {
	initmatdoub init_a = {{2, 1, -1, 1, -3}, {1, 0, 2, -1, 1}, {0, -2, -1, 1, -1}, {3, 1, -4, 0, 5}, {1, -1, -1, -1, 1}};
	initmatdoub init_b = {{7}, {2}, {-5}, {6}, {3}};
	matdoub a(5, 5, init_a);
	matdoub b(5, 1, init_b);
	gaussj(a, b);
	std::cout << b;
	return 0;
}