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
	initmatdoub init_a = {{4, -1, 1},
						  {2, 5, 2},
						  {1, 2, 4}};
	initvecdoub init_b = {8, 3, 11};
	matdoub a(3, 3, init_a);
	vecdoub b(3, init_b);
	LUdcmp lu(a);
	vecdoub x = lu.solve(b);
	std::cout << x;
	return 0;
}