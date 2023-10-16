import numpy as np


eps = 1e-3


def gen_m(n):
    return [[min(i+1, j+1)/max(i+1, j+1) for j in range(n)] for i in range(n)]


def write_m(m, filename):
    with open(filename, "w+") as f:
        f.writelines([f"{' '.join([str(e) for e in row])}\n" for row in m])


def read_m(filename):
    with open(filename) as f:
        return [map(int, row.strip('\n').split()) for row in f.readlines()]


def find_max_elem(S):
    max_elem = 0
    im, jm = 0, 0
    for i in range(S.shape[0] - 1):
        for j in range(i + 1, S.shape[0]):
            if abs(S[i, j]) > max_elem:
                max_elem, im, jm = abs(S[i, j]), i, j
    return im, jm


def rot_multiplication(A, cos, sin, i, j, side, herm):
    n = len(A)
    if side == 'l':
        ai = [A[i, k] * cos + A[j, k] * sin for k in range(n)]
        aj = [A[i, k] * (-sin) + A[j, k] * cos for k in range(n)]
        A[i] = ai
        A[j] = aj
    else:
        if herm:
            ai = [[A[k, i] * cos - A[k, j] * sin] for k in range(n)]
            aj = [[A[k, i] * sin + A[k, j] * cos] for k in range(n)]
        else:
            ai = [A[k, i] * cos - A[k, j] * sin for k in range(n)]
            aj = [A[k, i] * sin + A[k, j] * cos for k in range(n)]
        A[:, i] = ai
        A[:, j] = aj
    return A


def rotations(M, herm, tol=eps, max_iter=1000):
    A = np.matrix(M)
    S = A.copy()
    for i in range(max_iter):
        a, b = find_max_elem(S)
        if abs(S[a, b]) < tol:
            break
        theta = 0.5 * np.arctan2(2 * S[a, b], S[b, b] - S[a, a])
        c = np.cos(theta)
        s = np.sin(theta)
        S = rot_multiplication(rot_multiplication(S, c, -s, a, b, 'l', herm), c, s, a, b, 'r', herm)
    eigenvalues = np.diag(S)
    return eigenvalues


print(rotations(np.array([[1., 1., 1., 1.],
              [1., 2., 2., 2.],
              [1., 2., 3., 3.],
              [1., 2., 3., 4.]]), True))