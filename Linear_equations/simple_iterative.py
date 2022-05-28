import numpy as np

def jacobi(A, b, x0, tol, MAX_ITER):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x0 - n size column initial guess
    tol - tolerance value
    MAX_ITER - max iteration number
    x - n size column
    '''
    k = 1
    while k <= MAX_ITER:
        # formula suppose to sum A[i, j]*x[j] over j != i,
        # but there is a way to implement it without a loop,
        # if we modify A by making zeros on its diagonal
        A_mod = A - np.diag(np.diag(A))
        x = (-np.dot(A_mod, x0) + b)/np.diag(A)
        if np.linalg.norm(x - x0) < tol:
            return x
        k += 1
        x0 = x
    print('Maximum iterations exceeded')
    return

def gauss_seidel(A, b, x0, tol, MAX_ITER):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x0 - n size column initial guess
    tol - tolerance value
    MAX_ITER - max iteration number
    x - n size column
    '''
    k = 1
    n = b.shape[0]
    while k <= MAX_ITER:
        x = np.zeros(n)
        for i in range(n):
            x[i] = (-np.dot(A[i, :i], x[:i]) - np.dot(A[i, i+1:], x0[i+1:]) + b[i])/A[i, i]
        if np.linalg.norm(x - x0) < tol:
            return x
        k += 1
        x0 = x
    print('Maximum iterations exceeded')
    return

def SOR(A, b, x0, w, tol, MAX_ITER):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x0 - n size column initial guess
    w - is a SOR paramater (usually omega)
    tol - tolerance value
    MAX_ITER - max iteration number
    x - n size column
    '''
    k = 1
    n = b.shape[0]
    while k <= MAX_ITER:
        x = np.zeros(n)
        for i in range(n):
            x[i] = (1 - w)*x0[i] + w*(-np.dot(A[i, :i], x[:i]) - np.dot(A[i, i+1:], x0[i+1:]) + b[i])/A[i, i]
        if np.linalg.norm(x - x0) < tol:
            return x
        k += 1
        x0 = x
    print('Maximum iterations exceeded')
    return