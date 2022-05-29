import numpy as np

def gauss_elim(A, b):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    x - n size column
    '''
    n = b.shape[0]
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
        block[i+1:] -= np.outer(block[i+1:, i], block[i])/block[i, i]
    if block[-1, -2] == 0:
        print('No unique solution')
        return
    
    # backward substitution
    x = np.zeros(n)
    x[-1] = block[-1, -1]/block[-1, -2]
    for i in range(n-2, -1, -1):
        x[i] = (block[i, -1] - np.dot(block[i, i:-1], x[i:]))/block[i, i]

    return x