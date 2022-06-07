#ifndef GAUSSJ_H
#define GAUSSJ_H

#include "../myvector.hpp"
#include "../mymatrix.hpp"
#include "../util.hpp"

void gauss(matdoub &a, matdoub &b)
/*
Implementation of Gauss elimination with partial pivoting.
Parameters:
-----------
Input:
a: square (n,n) matrix
b: (n,m) matrix

Output:
b: (n,m) solution matrix

Runtime parameters:
irow: position of pivoting row
indxr: contains structure of initial a matrix, not necessary for solving
max: need to find maximum absolute value in column
pivinv: 1 / (pivot element)
temp: temporary variable
*/
{
    int irow, n = a.nrows(), m = b.ncols();
    double max, pivinv, temp;
    vecint indxr(n);
    for (int i = 0; i < n; ++i)
    {
        max = 0.0;
        for (int j = i; j < n; ++j)
        {
            if (fabs(a[j][i]) >= max)
            {
                max = fabs(a[j][i]);
                irow = j;
            }
        }

        if (irow != i)
        {
            for (int j = 0; j < n; ++j)
                SWAP(&a[irow][j], &a[i][j]);
            for (int j = 0; j < m; ++j)
                SWAP(&b[irow][j], &b[i][j]);
        }
        indxr[i] = irow;
        if (a[i][i] == 0.0)
            throw std::runtime_error("gauss: Singular matrix.");
        pivinv = 1.0 / a[i][i];
        for (int j = i + 1; j < n; ++j)
        {
            temp = a[j][i];
            for (int k = 0; k < n; ++k)
                a[j][k] -= a[i][k] * temp * pivinv;
            for (int k = 0; k < m; ++k)
                b[j][k] -= b[i][k] * temp * pivinv;
        }
    }

    for (int i = 0; i < m; ++i)
    {
        b[n - 1][i] /= a[n - 1][n - 1];
        for (int j = n - 2; j >= 0; --j)
        {
            for (int k = j + 1; k < n; ++k)
                b[j][i] -= a[j][k] * b[k][i];
            b[j][i] /= a[j][j];
        }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        if (indxr[i] != i)
            for (int j = 0; j < n; ++j)
                SWAP(&a[j][indxr[i]], &a[j][i]);
    }
}

void gauss(matdoub &a)
{
    matdoub b(a.nrows(), 1, 0);
    gauss(a, b);
}

#endif