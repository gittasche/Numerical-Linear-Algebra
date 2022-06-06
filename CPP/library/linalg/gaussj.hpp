#ifndef GAUSSJ_H
#define GAUSSJ_H

#include "../myvector.hpp"
#include "../mymatrix.hpp"
#include "../util.hpp"

void gaussj(matdoub &a, matdoub &b)
{
    int irow, icol, n = a.nrows(), m = b.ncols();
    vecint indxc(n), indxr(n), ipiv(n, 0);
    double max, pivinv, temp;
    for (int i = 0; i < n; ++i)
    {
        max = 0.0;
        for (int j = 0; j < n; ++j)
            if (ipiv[j] != 1)
                for (int k = 0; k < n; ++k)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j][k]) >= max)
                        {
                            max = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv[icol]);

        if (irow != icol)
        {
            for (int j = 0; j < n; ++j)
                SWAP(&a[irow][j], &a[icol][j]);
            for (int j = 0; j < m; ++j)
                SWAP(&b[irow][j], &b[icol][j]);
        }
        indxr[i] = irow;
        indxr[i] = icol;
        if (a[irow][icol] == 0.0)
            throw("gaussj: Singular matrix.");
        pivinv = 1.0/a[icol][irow];
        a[icol][irow] = 1.0;
        for (int j = 0; j < n; ++j)
            a[icol][j] *= pivinv;
        for (int j = 0; j < m; ++j)
            b[icol][j] *= pivinv;
        for (int j = 0; j < n; ++j)
            if (j != icol)
            {
                temp = a[j][icol];
                a[j][icol] = 0.0;
                for (int k = 0; k < n; ++k)
                    a[j][k] -= a[icol][k] * temp;
                for (int k = 0; k < m; ++k)
                    b[j][k] -= b[icol][k] * temp;
            }
    }

    for (int i = n - 1; i >= 0; --i)
    {
        if (indxr[i] != indxc[i])
            for (int j = 0; j < n; ++j)
                SWAP(&a[j][indxr[i]], &a[j][indxc[i]]);
    }
}

void gaussj(matdoub &a)
{
    matdoub b(a.nrows(), 1, 0);
    gaussj(a, b);
}

#endif