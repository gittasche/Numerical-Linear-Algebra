#include "../CPP/library/myvector.hpp"
#include "../CPP/library/mymatrix.hpp"
#include "../CPP/library/linalg/gaussj.hpp"

int main() {
	initmatdoub init_a = {{4, -1, 1}, {2, 5, 2}, {1, 2, 4}};
	initmatdoub init_b = {{8}, {3}, {11}};
	matdoub a(3, 3, init_a);
	matdoub b(3, 1, init_b);
	gaussj(a, b);
	std::cout << b;
	return 0;
}