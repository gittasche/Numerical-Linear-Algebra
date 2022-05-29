import numpy as np

# for inverse power
from _util import LU_solve

def power(A, x, tol, MAX_ITER):
    '''
    Power method for dominant eigenvalue
    and corresponding eigenvector finding
    Parameters:
    -----------
    A - (n, n) matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    k = 1
    p = np.argmax(np.abs(x))
    x /= x[p]
    while k <= MAX_ITER:
        y = np.dot(A, x)
        u = y[p]
        p = np.argmax(np.abs(y))
        if y[p] == 0.0:
            print ('Eigenvector x have eigenvalue 0')
            return x, 0
        err = np.linalg.norm(x - y/y[p], ord=np.inf)
        x = y/y[p]
        if err < tol:
            return x, u
        k += 1
    print('Maximum iterations exceeded')
    return

def aitken(A, x, tol, MAX_ITER):
    '''
    Aitken boosting for power method
    Parameters:
    -----------
    A - (n, n) matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    k = 1
    u0 = u1 = 0
    p = np.argmax(np.abs(x))
    x /= x[p]
    while k <= MAX_ITER:
        y = np.dot(A, x)
        u = y[p]
        u_ = u0 - (u1 - u0)**2/(u - u1 + u0)
        p = np.argmax(np.abs(y))
        if y[p] == 0.0:
            print ('Eigenvector x have eigenvalue 0')
            return x, 0
        err = np.linalg.norm(x - y/y[p], ord=np.inf)
        x = y/y[p]
        if err < tol and k >= 4:
            return x, u_
        k += 1
        u0 = u1
        u1 = u
    print('Maximum iterations exceeded')
    return

def symm_power(A, x, tol, MAX_ITER):
    '''
    Symmetric power method for dominant eigenvalue
    and corresponding eigenvector finding
    Parameters:
    -----------
    A - (n, n) symmetric matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    k = 1
    x /= np.linalg.norm(x)
    while k <= MAX_ITER:
        y = np.dot(A, x)
        u = np.dot(x, y)
        if np.linalg.norm(y) == 0.0:
            print ('Eigenvector x have eigenvalue 0')
            return x, 0
        err = np.linalg.norm(x - y/np.linalg.norm(y))
        x = y/np.linalg.norm(y)
        if err <= tol:
            return x, u
        k += 1
    print('Maximum iterations exceeded')
    return

def aitken_symm(A, x, tol, MAX_ITER):
    '''
    Aitken boosting for symmetric power method
    Parameters:
    -----------
    A - (n, n) symmetric matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    k = 1
    u0 = u1 = 0
    x /= np.linalg.norm(x)
    while k <= MAX_ITER:
        y = np.dot(A, x)
        u = np.dot(x, y)
        u_ = u0 - (u1 - u0)**2/(u - u1 + u0)
        if np.linalg.norm(y) == 0.0:
            print ('Eigenvector x have eigenvalue 0')
            return x, 0
        err = np.linalg.norm(x - y/np.linalg.norm(y))
        x = y/np.linalg.norm(y)
        if err <= tol and k >= 4:
            return x, u_
        k += 1
        u0 = u1
        u1 = u
    print('Maximum iterations exceeded')
    return

def inv_power(A, x, tol, MAX_ITER):
    '''
    Inverse power method for dominant eigenvalue
    and corresponding eigenvector finding
    Parameters:
    -----------
    A - (n, n) matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    n = A.shape[0]
    q = x.T @ A @ x / np.dot(x, x)
    k = 1
    p = np.argmax(np.abs(x))
    x /= x[p]
    while k <= MAX_ITER:
        y = LU_solve(A - q*np.identity(n), x)
        if y is None:
            return None, q
        u = y[p]
        p = np.argmax(np.abs(y))
        err = np.linalg.norm(x - y/y[p], ord=np.inf)
        x = y/y[p]
        if err < tol:
            u = 1/u + q
            return x, u
        k += 1
    print('Maximum iterations exceeded')
    return

def aitken_inv(A, x, tol, MAX_ITER):
    '''
    Aitken boosting for inverse power method
    Parameters:
    -----------
    A - (n, n) matrix
    x - n size column approximate eigenvector
    tol - tolerance value
    MAX_ITER - maximum iterations number
    '''
    n = A.shape[0]
    q = x.T @ A @ x / np.dot(x, x)
    k = 1
    u0 = u1 = 0
    p = np.argmax(np.abs(x))
    x /= x[p]
    while k <= MAX_ITER:
        y = LU_solve(A - q*np.identity(n), x)
        if y is None:
            return None, q
        u = y[p]
        u_ = u0 - (u1 - u0)**2/(u - u1 + u0)
        p = np.argmax(np.abs(y))
        err = np.linalg.norm(x - y/y[p], ord=np.inf)
        x = y/y[p]
        if err < tol and k >= 4:
            u_ = 1/u_ + q
            return x, u_
        k += 1
        u0 = u1
        u1 = u
    print('Maximum iterations exceeded')
    return