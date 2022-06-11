#include "../CPP/library/myvector.hpp"
#include "../CPP/library/mymatrix.hpp"
#include "../CPP/library/linalg/gauss.hpp"
#include "../CPP/library/linalg/eigen_sym.hpp"
#include "../CPP/library/linalg/ludcmp.hpp"

int main()
{
	initmatdoub init_a = {{4, -1, 1},
						  {-1, 3, -2},
						  {1, -2, 3}};
	matdoub a(3, 3, init_a);
	Jacobi jac(a, 50);
	jac.solve();
	jac.print_res('Y');
	return 0;
}