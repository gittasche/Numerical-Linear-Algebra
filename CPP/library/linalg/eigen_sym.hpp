#ifndef EIGSYM_H
#define EIGSYM_H

#include "../myvector.hpp"
#include "../mymatrix.hpp"
#include "../util.hpp"

void rotate(matdoub &mat, const double s, const double tau, const int i, const int j, const int k, const int l);

class Jacobi
{
private:
    const int n;
    matdoub a, vecs;
    vecdoub lams;
    int nrot;
    const double EPS;
    const int MAX_ITER;
    bool solved = false;

public:
    Jacobi(matdoub &aa, const int MAX_ITER);
    void solve();
    friend void rotate(matdoub &mat, const double s, const double tau, const int i, const int j, const int k, const int l);
    void eigsrt(const char vecsort = 'Y');
    void print_res(const char vecprint = 'N');
};

Jacobi::Jacobi(matdoub &aa, const int MAX_ITER) : n(aa.nrows()), a(aa), vecs(n, n), lams(n), nrot(0), EPS(std::numeric_limits<double>::epsilon()), MAX_ITER(MAX_ITER) {}

void Jacobi::solve()
/*
Implementation of Jacobi rotation method
for finding eigenvectors and eigenvalues of
symmetric matrix.

Parameters:
-----------
Input:
a: symmetric square (n,n) matrix
MAX_ITER: number of maximum iterations

Output:
vecs: matrix of normalized eigenvectors in columns
lams: vector of eigenvectors

Runtime parameters:
sum: sum of above diagonal elements
tresh: threshold value for convergence boosting
g: off-diagonal element to check its magnitude
h: temporary variable
t: root of equation t^2 + 2*t*theta - 1 = 0
theta: cot(2 * rotation_angle)
c: diagonal coefs of rotation matrix (cos(rotation_angle))
s: off-diagonal coefs of rotation matrix (sin(rotation_angle))
tau: tan(rotation_angle / 2) equals s / (1 + c)
b: temp vector
z: temp vector
*/
{
    double sum, tresh, g, h, t, theta, c, s, tau;
    vecdoub b(n), z(n);
    vecs.identity();
    for (int ip = 0; ip < n; ++ip)
    {
        b[ip] = lams[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    for (int i = 1; i <= MAX_ITER; ++i)
    {
        sum = 0.0;
        for (int ip = 0; ip < n - 1; ++ip)
            for (int iq = ip + 1; iq < n; ++iq)
                sum += fabs(a[ip][iq]);
        if (sum == 0.0)
        {
            this->eigsrt();
            solved = true;
            return;
        }
        if (i < 4)
            tresh = 0.2 * sum / (n * n);
        else
            tresh = 0.0;
        for (int ip = 0; ip < n - 1; ++ip)
        {
            for (int iq = ip + 1; iq < n; ++iq)
            {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && g <= EPS * fabs(lams[ip]) && g <= EPS * fabs(lams[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh)
                {
                    h = lams[iq] - lams[ip];
                    if (g <= EPS * fabs(h))
                        t = a[ip][iq] / h;
                    else
                    {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + SQR(theta)));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + SQR(t));
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    lams[ip] -= h;
                    lams[iq] += h;
                    a[ip][iq] = 0.0;
                    for (int j = 0; j < ip; ++j)
                        rotate(a, s, tau, j, ip, j, iq);
                    for (int j = ip + 1; j < iq; ++j)
                        rotate(a, s, tau, ip, j, j, iq);
                    for (int j = iq + 1; j < n; ++j)
                        rotate(a, s, tau, ip, j, iq, j);
                    for (int j = 0; j < n; ++j)
                        rotate(vecs, s, tau, j, ip, j, iq);
                    ++nrot;
                }
            }
        }
        for (int ip = 0; ip < n; ++ip)
        {
            b[ip] += z[ip];
            lams[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    std::cout << "Maximum iterations exceeded." << std::endl;
}

void rotate(matdoub &mat, const double s, const double tau, const int i, const int j, const int k, const int l)
{
    double g = mat[i][j];
    double h = mat[k][l];
    mat[i][j] = g - s * (h + g * tau);
    mat[k][l] = h + s * (g - h * tau);
}

void Jacobi::eigsrt(const char vecsort)
/*
Sorting eigenvalues and corresponding eigenvectors
in descending order (insertion algorithm).

Parameters:
-----------
vecsort: 'Y' if you need to sort eigenvectors
*/
{
    int k, n = lams.size();
    double p;
    for (int i = 0; i < n - 1; ++i)
    {
        k = i;
        p = lams[k];
        for (int j = i; j < n; ++j)
            if (lams[j] >= p)
            {
                k = j;
                p = lams[k];
            }
        if (k != i)
        {
            lams[k] = lams[i];
            lams[i] = p;
            if (vecsort == 'Y')
            {
                for (int j = 0; j < n; ++j)
                    SWAP(&vecs[j][i], &vecs[j][k]);
            }
        }
    }
}

void Jacobi::print_res(const char vecprint)
/*
Print results

Parameters:
-----------
vecprint: 'Y' if you need to print eigenvectors
*/
{
    if (solved)
    {
        std::cout << "Eigenvalues:\n"
                  << lams << "\n";
        if (vecprint == 'Y')
            std::cout << "Normalized eigenvectors in columns:\n"
                      << vecs << std::endl;
    }
    else
        std::cout << "Not solved." << std::endl;
}
#endif