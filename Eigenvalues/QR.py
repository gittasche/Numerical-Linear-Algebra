import numpy as np

def QR_decomp(a, b, tol, MAX_ITER):
    #TODO: impement calculation of Q matrix
    '''
    QR method for finding the eigenvalues
    of (n,n) symmetric tridiagonal matrix.
    This code doesn't return a Q matrix.
    Use Householder method to tridiagonalize
    matrix if ist not already.

    Parameters:
    -----------
    a - main diagonal of input matrix
    b - adjecent diagonals to a
    tol - tolerance value
    MAX_ITER - maximum iteration number

    If you have tridiagonal matrix A put
    a = np.diag(A) and b = np.diag(A, k=1)
    '''
    a_k = np.copy(a)
    b_k = np.copy(b)
    n = a_k.shape[0]
    k = 1
    SHIFT = 0
    lams = np.array([])


    while k <= MAX_ITER:

        if np.abs(b_k[-1]) < tol:
            lams = np.append(lams, a_k[-1] + SHIFT)
            a_k = np.delete(a_k, -1)
            b_k = np.delete(b_k, -1)
            n -= 1
            if n == 0:
                return lams, R
        
        if np.abs(b_k[0]) < tol:
            lams = np.append(lams, a_k[0] + SHIFT)
            a_k = np.delete(a_k, 0)
            b_k = np.delete(b_k, 0)
            n -= 1
            if n == 0:
                return lams, R

        if n == 1:
            lams = np.append(lams, a_k[0] + SHIFT)
            return lams, R

        b = -(a_k[-2] + a_k[-1])
        c = a_k[-1]*a_k[-2] - b_k[-1]**2
        d = np.sqrt(b**2 - 4*c)
        if b > 0:
            u = np.array([-2*c/(b + d), -(b + d)/2])
        else:
            u = np.array([(d - b)/2, 2*c/(d - b)])

        if n == 2:
            lams = np.append(lams, u + SHIFT)
            return lams, R

        sigma = u[np.argmin(np.abs(u - a_k[-1]))]
        SHIFT += sigma
        d = a_k - sigma


        # temp arrays
        x = np.zeros(n)
        y = np.zeros(n)

        # elements of rotation matrices
        g = np.zeros(n)
        s = np.zeros(n)

        # diags of R matrix
        z = np.zeros(n)
        q = np.zeros(n-1)
        r = np.zeros(n-2)

        x[0] = d[0]
        y[0] = b_k[0]
        for j in range(1, n):
            z[j-1] = np.sqrt(x[j-1]**2 + b_k[j-1]**2)
            g[j] = x[j-1]/z[j-1]
            s[j] = b_k[j-1]/z[j-1]
            q[j-1] = g[j]*y[j-1] + s[j]*d[j]
            x[j] = -s[j]*y[j-1] + g[j]*d[j]
            if j != n-1:
                r[j-1] = s[j]*b_k[j]
                y[j] = g[j]*b_k[j]
        z[-1] = x[-1]

        # modifying main diagonal a
        a_k[0] = s[1]*q[0] + g[1]*z[0]
        a_k[1:-1] = s[2:]*q[1:] + g[1:-1]*g[2:]*z[1:-1]
        a_k[-1] = g[-1]*z[-1]

        # modifying adjecent diagonals b
        b_k = s[1:]*z[1:]

        # R matrix to return
        R = np.diag(z) + np.diag(q, k=1) + np.diag(r, k=2)

        k += 1
    print('Maximum iterations exceeded')
    return

def QR_parts(a, b, tol, MAX_ITER):
    '''
    QR method for eigenvalues, but it returns
    generator to obtain eigenvalues by parts.
    (Purpose of this function is a yield practice)

    QR method for finding the eigenvalues
    of (n,n) symmetric tridiagonal matrix.
    This code doesn't return a Q and R matrices.
    Use Householder method to tridiagonalize
    matrix if ist not already.

    Parameters:
    -----------
    a - main diagonal of input matrix
    b - adjecent diagonals to a
    tol - tolerance value
    MAX_ITER - maximum iteration number

    If you have tridiagonal matrix A put
    a = np.diag(A) and b = np.diag(A, k=1)
    '''
    a_k = np.copy(a)
    b_k = np.copy(b)
    n = a_k.shape[0]
    k = 1
    SHIFT = 0


    while k <= MAX_ITER:

        if np.abs(b_k[-1]) < tol:
            lam = a_k[-1] + SHIFT
            a_k = np.delete(a_k, -1)
            b_k = np.delete(b_k, -1)
            yield lam
            n -= 1
            if n == 0:
                return
        
        if np.abs(b_k[0]) < tol:
            lam = a_k[0] + SHIFT
            a_k = np.delete(a_k, 0)
            b_k = np.delete(b_k, 0)
            yield lam
            n -= 1
            if n == 0:
                return

        if n == 1:
            lam = a_k[0] + SHIFT
            yield lam
            return

        b = -(a_k[-2] + a_k[-1])
        c = a_k[-1]*a_k[-2] - b_k[-1]**2
        d = np.sqrt(b**2 - 4*c)
        if b > 0:
            u = np.array([-2*c/(b + d), -(b + d)/2])
        else:
            u = np.array([(d - b)/2, 2*c/(d - b)])

        if n == 2:
            lams = u + SHIFT
            yield lams[0]
            yield lams[1]
            return

        sigma = u[np.argmin(np.abs(u - a_k[-1]))]
        SHIFT += sigma
        d = a_k - sigma


        # temp arrays
        x = np.zeros(n)
        y = np.zeros(n)

        # elements of rotation matrices
        g = np.zeros(n)
        s = np.zeros(n)

        # diags of R matrix
        z = np.zeros(n)
        q = np.zeros(n-1)
        r = np.zeros(n-2)

        x[0] = d[0]
        y[0] = b_k[0]
        for j in range(1, n):
            z[j-1] = np.sqrt(x[j-1]**2 + b_k[j-1]**2)
            g[j] = x[j-1]/z[j-1]
            s[j] = b_k[j-1]/z[j-1]
            q[j-1] = g[j]*y[j-1] + s[j]*d[j]
            x[j] = -s[j]*y[j-1] + g[j]*d[j]
            if j != n-1:
                r[j-1] = s[j]*b_k[j]
                y[j] = g[j]*b_k[j]
        z[-1] = x[-1]

        # modifying main diagonal a
        a_k[0] = s[1]*q[0] + g[1]*z[0]
        a_k[1:-1] = s[2:]*q[1:] + g[1:-1]*g[2:]*z[1:-1]
        a_k[-1] = g[-1]*z[-1]

        # modifying adjecent diagonals b
        b_k = s[1:]*z[1:]

        k += 1
    print('Maximum iterations exceeded')
    return