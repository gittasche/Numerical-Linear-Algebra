import numpy as np

def conj_grads(A, b, C, x0, tol, MAX_ITER):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    C - preconditioning matrix (usually identity)
    x0 - initial guess n size column
    tol - tolerance value
    MAX_ITER - max iteration number
    x - n size column
    '''
    x = x0
    C_inv = np.linalg.inv(C)
    r = b - np.dot(A, x)
    w = np.dot(C_inv, r)
    v = np.dot(C_inv.T, w)
    alpha = np.dot(w, w)

    k = 1
    while k <= MAX_ITER:
        if np.linalg.norm(v) < tol:
            return x, r

        u = np.dot(A, v)
        t = alpha/np.dot(u, v)
        x += t*v
        r -= t*u
        w = np.dot(C_inv, r)
        beta = np.dot(w, w)

        if beta < tol and np.linalg.norm(r) < tol:
            return x, r

        s = beta/alpha
        v = np.dot(C_inv.T, w) + s*v
        alpha = beta
        k += 1

    print('Maximum iterations exceeded')
    return