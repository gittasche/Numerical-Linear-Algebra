#ifndef LUDCMP_H
#define LUDCMP_H

#include "../myvector.hpp"
#include "../mymatrix.hpp"
#include "../util.hpp"

class LUdcmp
{
private:
    int n;
    matdoub lu;
    vecint indx;
    double d;

public:
    LUdcmp(matdoub &a);
    vecdoub solve(vecdoub &b);
    matdoub solve(matdoub &b);
    matdoub inverse();
    double det();
};

LUdcmp::LUdcmp(matdoub &a) : n(a.nrows()), lu(a), indx(n)
/*
Implementation of LU decomposition.
Parameters:
-----------
Input:
a: square (n,n) matrix to be decomposed

Output:
lu: square (n,n) matrix contains both L and U matrices (Crout method)

Runtime parameters:
TINY: tiny double number which will be denominator instead of zero
jmax: index of pivoting element
max: pivoting element
vv: temp vector of 1.0/max's
*/
{
    const double TINY = 1.0e-40;
    int jmax;
    double max;
    vecdoub vv(n);
    d = 1.0;

    for (int i = 0; i < n; ++i)
    {
        max = 0.0;
        for (int j = 0; j < n; ++j)
            if (fabs(lu[i][j]) > max)
                max = fabs(lu[i][j]);
        if (max == 0.0)
            throw std::runtime_error("LUdcmp: Singular matrix.");
        vv[i] = 1.0 / max;
    }

    for (int i = 0; i < n; ++i)
    {
        max = 0.0;
        for (int j = i; j < n; ++j)
            if (vv[i] * fabs(lu[j][i]) > max)
            {
                max = vv[i] * fabs(lu[j][i]);
                jmax = i;
            }
        if (i != jmax)
        {
            for (int j = 0; j < n; ++j)
                SWAP(&lu[jmax][j], &lu[i][j]);
            d = -d;
            vv[jmax] = vv[i];
        }

        indx[i] = jmax;
        if (lu[i][i] == 0.0)
            lu[i][i] = TINY;

        for (int j = i + 1; j < n; ++j)
        {
            lu[j][i] /= lu[i][i];
            for (int k = i + 1; k < n; ++k)
                lu[j][k] -= lu[j][i] * lu[i][k];
        }
    }
    std::cout << lu << std::endl;
}

vecdoub LUdcmp::solve(vecdoub &b)
{
    if (n != b.size())
        throw std::runtime_error("LUdcmp: Non-equal sizes.");

    int ip, ii = 0;
    double sum;
    vecdoub x(b);
    for (int i = 0; i < n; ++i)
    {
        ip = indx[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (int j = ii - 1; j < i; ++j)
                sum -= lu[i][j] * x[j];
        else if (sum != 0.0)
            ii = i + 1;
        x[i] = sum;
    }

    for (int i = n - 1; i >= 0; --i)
    {
        sum = x[i];
        for (int j = i + 1; j < n; ++j)
            sum -= lu[i][j] * x[j];
        x[i] = sum / lu[i][i];
    }

    return x;
}

matdoub LUdcmp::solve(matdoub &b)
{
    int m = b.ncols();
    if (n != b.nrows())
        throw std::runtime_error("LUdcmp: Non-equal sizes.");
    
    vecdoub xx(n);
    matdoub x(n, m);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
            xx[j] = b[j][i];
        xx = LUdcmp::solve(xx);
        for (int j = 0; j < n; ++j)
            x[j][i] = xx[j];
    }

    return x;
}

matdoub LUdcmp::inverse()
{
    matdoub ainv(n, n);
    ainv.identity();
    ainv = solve(ainv);
    return ainv;
}

double LUdcmp::det()
{
    double dd = d;
    for (int i = 0; i < n; ++i)
        dd *= lu[i][i];
    return dd;
}

#endif