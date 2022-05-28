import numpy as np
# iterative refinement requires Gauss elimination
from gauss_elim import gauss_elim

def refinement(A, b, tol, MAX_ITER, t):
    '''
    A is (n, n) non-singular matrix
    Equation Ax = b
    Parameters:
    ------------
    A - square (n, n) size matrix
    b - n size column
    tol - tolerance value
    MAX_ITER - max iteration number
    t - digit of precision
    x - n size column
    '''
    k = 1
    x = gauss_elim(A, b)
    while k <= MAX_ITER:
        r = b - np.dot(A, x)
        y = gauss_elim(A, r)
        xx = x + y
        if k == 1:
            COND = np.linalg.norm(y, ord=np.inf)/np.linalg.norm(xx, ord=np.inf) * 10**t
        if np.linalg.norm(xx - x, ord=np.inf) < tol:
            return xx, COND
        k += 1
    print('Maximum iterations exceeded')
    return COND
def main():
    # these arrays have low precision float16 type
    A = np.array([[3.3330, 15920, -10.333],
                [2.2220, 16.710, 9.6120],
                [1.5611, 5.1791, 1.6852]], dtype=np.float16)
    b = np.array([15913, 28.544, 8.4254], dtype=np.float16)

    x, COND = refinement(A, b, 1e-3, 200, 5)
    print(f'Solution: {x}\nCondition number: {COND}')

if __name__ == '__main__':
    main()