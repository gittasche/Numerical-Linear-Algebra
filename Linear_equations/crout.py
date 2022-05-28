import numpy as np

def crout(A, b):
    '''
    A is tridiagonal (n, n) matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x - n size column
    '''
    n = b.shape[0]
    block = np.c_[A, b].astype(np.float64)
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    z = np.zeros(n)

    L[0, 0] = block[0, 0]
    U[0, 1] = block[0, 1]/L[0, 0]
    z[0] = block[0, -1]/L[0, 0]

    for i in range(1, n-1):
        L[i, i-1] = block[i, i-1]
        L[i, i] = block[i, i] - L[i, i-1]*U[i-1, i]
        U[i, i+1] = block[i, i+1]/L[i, i]
        z[i] = (block[i, -1] - L[i, i-1]*z[i-1])/L[i, i]

    L[-1, -2] = block[-1, -2]
    L[-1, -1] = block[-1, -1] - L[-1, -2]*U[-2, -1]
    z[-1] = (block[-1, -1] - L[-1, -2]*z[-2])/L[-1, -1]

    # solve Ux = z
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = z[i] - U[i, i+1]*x[i+1]

    return x