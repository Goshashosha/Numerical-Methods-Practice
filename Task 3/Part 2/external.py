import numpy as np


def gen_m(n):
    return np.array([[min(i+1, j+1) for j in range(n)] for i in range(n)], float)


def write_m(m, filename):
    with open(filename, "w+") as f:
        f.writelines([f"{' '.join(row)}\n" for row in m])


def read_m(filename):
    with open(filename) as f:
        return [map(int, row.strip('\n').split()) for row in f.readlines()]


def hh_trans(a):
    n = a.shape[0]
    sx = np.zeros(n)

    for i in range(n - 2):
        for k in range(i, n):
            sx[k] = a[n - 1, i] * a[n - 1, k]

        for j in range(n - 2, i, -1):
            sx[i] = sx[i] + a[j, i] * a[j, i]

        for k in range(i + 1, n):
            for j in range(n - 2, k - 1, -1):
                sx[k] = sx[k] + a[j, i] * a[j, k]

            for j in range(k - 1, i, -1):
                sx[k] = sx[k] + a[j, i] * a[k, j]

        alpha = np.sqrt(sx[i])

        if a[i + 1, i] != 0:
            beta = 1. / alpha

            for j in range(i + 2, n):
                a[j, i] = a[j, i] * beta

            sx[i] = a[i + 1, i] * beta + np.sign(a[i + 1, i])
            a[i + 1, i] = alpha
            g = 1 / abs(sx[i])
            sx2 = 0

            for k in range(i + 2, n):
                sx[k] = sx[k] * beta * g + np.sign(sx[i]) * abs(a[k, i + 1])
                sx2 = sx[k] * a[k, i] + sx2

            sx2 = g * sx2

            for k in range(i + 2, n):
                a[k, k] = a[k, k] - 2 * a[k, i] * sx[k] + sx2 * a[k, i] ** 2

                for j in range(k + 1, n):
                    a[j, k] = a[j, k] - a[j, i] * sx[k] - a[k, i] * sx[i] + sx2 * a[j, i] * a[k, i]
        else:
            if alpha != 0:
                beta = 1 / alpha

                for j in range(i + 2, n):
                    a[j, i] = a[j, i] * beta

                sx[i] = -1.
                a[i + 1, i] = alpha
                g = 1.
                sx2 = 0.

                for k in range(i + 2, n):
                    sx[k] = sx[k] * beta * g + np.sign(sx[i]) * abs(a[k, i + 1])
                    sx2 = sx[k] * a[k, i] + sx2

                sx2 = g * sx2

                for k in range(i + 2, n):
                    a[k, k] = a[k, k] - 2 * a[k, i] * sx[k] + sx2 * a[k, i] ** 2

                    for j in range(k + 1, n):
                        a[j, k] = a[j, k] - a[j, i] * sx[k] - a[k, i] * sx[i] + sx2 * a[j, i] * a[k, i]
            else:
                sx[i] = 1.

        return a.diagonal(), a.diagonal(-1)


