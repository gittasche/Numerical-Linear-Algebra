import numpy as np

def gauss_elim(A, b):
    n = len(b)
    block = np.c_[A, b].astype(np.float64)
    for i in range(n - 1):
        p = np.argmax(np.abs(block[i:, i])) + i
        if block[p, i] == 0:
            print('No unique solution')
            return
        
        if p != i:
            block[[i, p]] = block[[p, i]]
        
        block[i+1:] -= block[i+1:, i, np.newaxis] @ block[np.newaxis, i]/block[i, i]
        
    if block[-1, -2] == 0:
        print('No unique solution')
        return
        
    x = np.zeros(n)
    x[-1] = block[-1, -1]/block[-1, -2]
    for i in range(n-2, -1, -1):
        x[i] = (block[i, -1] - np.dot(block[i, i:-1], x[i:]))/block[i, i]

    return x