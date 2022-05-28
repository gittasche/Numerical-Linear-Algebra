import numpy as np

def doolittle(A):
    '''
    Decomposition A = LU,
    where L left-bottom triange with ones on diagonal
    and U right-upper triangle
    Parameters:
    ------------
    A - square (n, n) size matrix
    '''
    n = A.shape[0]
    L = np.identity(n)
    U = np.zeros((n, n))

    if A[0, 0] == 0:
        print('Factorization impossible')
        return
    
    U[0, :] = A[0, :]/L[0, 0]
    L[1:, 0] = A[1:, 0]/U[0, 0]

    for i in range(1, n):
        U[i, i:] = A[i, i:] - np.dot(L[i, :i], U[:i, i:])
        L[i+1:, i] = (A[i+1:, i] - np.dot(L[i+1:, :i], U[:i, i]))/U[i, i]
        if U[i, i] == 0:
            print('Factorization impossible')
            return

    return L, U

def LU_solve(A, b):
    '''
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x - n size column
    '''
    n = b.shape[0]
    L, U = doolittle(A)

    # forward substitution Ly = b,
    # where y = Ux
    y = np.zeros(n)
    y[0] = b[0]/L[0, 0]
    for i in range(1, n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i]))/L[i, i]
    
    # backward substitution Ux = y
    x = np.zeros(n)
    x[-1] = y[-1]/U[-1, -1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i:], x[i:]))/U[i, i]
    
    return x