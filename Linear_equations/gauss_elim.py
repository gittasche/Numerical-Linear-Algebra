import numpy as np

def gauss_elim(A, b):
    n = len(b)
    block = np.c_[A, b]
    for i in range(n - 1):
        p = np.argmax(np.abs(A[i:, i])) + i
        if block[p, i] == 0:
            print('No unique solution')
            return
        
        if p != i:
            A[[i, p]] = A[[p, i]]
        
        for j in range(i+1, n):
            block[j] -= block[i]*block[j, i]/block[i, i]
        
    if block[-1, -2] == 0:
        print('No unique solution')
        return
        
    x = np.zeros(n)
    x[-1] = block[-1, -1]/block[-1, -2]
    for i in range(n-2, -1, -1):
        x[i] = (block[i, -1] - np.dot(block[i, i:-1], x[i:]))/block[i, i]

    return x
