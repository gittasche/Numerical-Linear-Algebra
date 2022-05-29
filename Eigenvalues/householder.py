import numpy as np

def householder(A):
    '''
    Householder algorithm to transform
    square (n,n) matrix into tridiagonal
    form.
    This method doesn't solve eigenproblem
    by itself. It is useful for QR decomposition
    Parameters:
    -----------
    A - (n,n) size matrix
    '''
    n = A.shape[0]
    v = np.zeros(n)
    u = np.zeros(n)
    z = np.zeros(n)

    for k in range(n-2):
        q = np.sum(A[k+1:, k]**2)
        if A[k+1, k] == 0:
            alpha = -np.sqrt(q)
        else:
            alpha = -np.sqrt(q)*np.sign(A[k+1, k])
        RSQ = alpha**2 - alpha*A[k+1, k]
        v[k] = 0
        v[k+1] = A[k+1, k] - alpha
        v[k+2:] = A[k+2:, k]
        u[k:] = np.dot(A[k:, k+1:], v[k+1:])/RSQ
        PROD = np.dot(v[k+1:], u[k+1:])
        z[k:] = u[k:] - (PROD/2/RSQ)*v[k:]
        A -= np.outer(v, z) + np.outer(z, v)
    return A