import numpy as np

def LDLt(A):
    '''
    A is positive define (n, n) matrix
    Decomposition A = LDLt,
    where L left-bottom triange with ones on diagonal,
    Lt is transposed L
    and D is diagonal
    Parameters:
    ------------
    A - square (n, n) size matrix
    '''
    n = A.shape[0]
    L = np.identity(n)
    D = np.zeros(n)

    for i in range(n):
        V = L[i, :i]*D[:i]
        D[i] = A[i, i] - np.dot(L[i, :i], V)
        L[i+1:, i] = (A[i+1:, i] - np.dot(L[i+1:, :i], V))/D[i]
    return L, D

def LDLt_solve(A, b):
    '''
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x - n size column
    '''
    n = b.shape[0]
    L, D = LDLt(A)

    # forward substitution Ly = b,
    # where y = DLtx
    z = np.zeros(n)
    z[0] = b[0]/L[0, 0]
    for i in range(1, n):
        z[i] = (b[i] - np.dot(L[i, :i], z[:i]))/L[i, i]

    # solve Dz = y, where z = Ltx
    y = z/D

    # backward substitution Ltx = z
    U = L.T
    x = np.zeros(n)
    x[-1] = y[-1]/U[-1, -1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i:], x[i:]))/U[i, i]

    return x
