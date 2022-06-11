#include "../CPP/library/myvector.hpp"
#include "../CPP/library/mymatrix.hpp"
#include "../CPP/library/linalg/gauss.hpp"
#include "../CPP/library/linalg/eigen_sym.hpp"
#include "../CPP/library/linalg/ludcmp.hpp"

int main()
{
	// Jacobi test
	// initmatdoub init_a = {{4, -1, 1},
	// 					  {-1, 3, -2},
	// 					  {1, -2, 3}};
	// matdoub a(3, 3, init_a);
	// Jacobi jac(a, 50);
	// jac.solve();
	// jac.print_res('Y');

	// LU test
	initmatdoub init_a = {{2, 0, 0, 0},
						  {1, 1.5, 0, 0},
						  {0, -3, 0.5, 0},
						  {2, -2, 1, 1}};
	matdoub a(4, 4, init_a);
	initvecdoub init_b = {3, 4.5, -6.6, 0.8};
	vecdoub b(4, init_b);
	LUdcmp lu(a);
	vecdoub x = lu.solve(b);
	std::cout << x;
	return 0;
}