import numpy as np

def gauss_jordan(A, b):
    '''
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x - n size column
    '''
    n = len(b)
    block = np.c_[A, b].astype(np.float64)
    for i in range(n - 1):

        # first element need to be nonzero
        p = np.argmax(np.abs(block[i:, i])) + i
        if block[p, i] == 0:
            print('No unique solution')
            return
        
        if p != i:
            block[[i, p]] = block[[p, i]]
        
        # elimination process
        block[i+1:] -= block[i+1:, i, np.newaxis] @ block[np.newaxis, i]/block[i, i]
        
    if block[-1, -2] == 0:
        print('No unique solution')
        return

    # diagonalizing matrix
    for i in range(n-2, -1, -1):
        coefs = np.array([block[i, k]/block[k, k] for k in range(i+1, n)])
        block[i, -1] -= np.dot(coefs, block[i+1:, -1])
    
    block[:, :-1] = np.diag(np.diag(block[:, :-1]))
    x = block[:, -1]/np.diag(block[:, :-1])

    return x