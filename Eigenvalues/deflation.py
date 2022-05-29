import numpy as np

from _util import aitken

def wielandt(A, lam, v, x, tol, MAX_ITER):
    '''
    Wielandt deflation method for second dominant
    eigenvalue and corresponding eigenvector finding
    Parameters:
    -----------
    A - (n, n) matrix
    lam - A dominant eigenvector
    v - lam eigenvector
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    i = np.argmax(np.abs(v))

    temp = A[i]/v[i]/lam
    B = A - lam*np.outer(v, temp)
    B = np.delete(np.delete(B, i, axis=1), i, axis=0)

    w, u = aitken(B, x, tol, MAX_ITER)
    if w is None:
        print('Method failed')
        return

    w = np.concatenate([w[:i], [0], w[i:]])
    w = (u - lam)*w + np.dot(A[i], w)*v/v[i]

    return w, u