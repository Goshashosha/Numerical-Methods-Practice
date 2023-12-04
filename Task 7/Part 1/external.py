import pandas as pd
import numpy as np


trap_ints = {'8': [], '16': [], '32': [], '64': []}


def k(x, s):
    return abs(x - s)


def f(x):
    return x ** 3


# ok
def linspace(a, b, n):
    return [a + i * (b - a) / (n - 1) for i in range(n)]


def zero_mat(n):
    return [[0. for _ in range(n)] for _ in range(n)]


# lambda
lbd = 1


def phi(i, x, linsps):
    if i == 0:
        return 1 - (x - linsps[i+1]) / (linsps[i+1] - linsps[i])
    elif i == len(linsps) - 1:
        return (x - linsps[i-1]) / (linsps[i] - linsps[i-1])
    else:
        return (x - linsps[i-1]) / (linsps[i] - linsps[i-1]) - (x - linsps[i+1]) / (linsps[i+1] - linsps[i])


# ok
def trap_int(x, func):
    s = sum([(x[i + 1] - x[i]) / 2 * (func[i] + func[i + 1]) for i in range(len(x) - 1)])
    trap_ints[f'{len(x)}'].append(s)
    return s


def trap_int_2(x, func):
    return sum([(x[i + 1] - x[i]) * func(x[i]) for i in range(len(x) - 1)])


def double_int(x, s, func):
    res = 0
    delta_x = x[1] - x[0]
    delta_s = s[1] - s[0]
    for i in range(len(x) - 1):
        for j in range(len(s) - 1):
            res += func[i][j] * delta_s * delta_x
    return res


def alhpa_mat(phi, x, n):
    a = zero_mat(n)
    for i in range(n):
        for j in range(n):
            a[i][j] = trap_int(x, [phi(i, point, x) * phi(j, point, x) for point in x])
    return a


def beta_mat(phi, k, x, n):
    b = zero_mat(n)
    for i in range(n):
        for j in range(n):
            b[i][j] = double_int(x, x, [[k(xx, ss) * phi(i, xx, x) * phi(j, xx, x) for ss in x] for xx in x])
    return b


def gamma_vec(phi, k, f, x, n):
    g = [0. for _ in range(n)]
    for i in range(n):
        g[i] = double_int(x, x, [[k(xx, ss) * phi(i, xx, x) * f(ss) for ss in x] for xx in x])
    return g


def matrix_sub(a, b):
    return [[a[i][j] - b[i][j] for i in range(len(a))] for j in range(len(a))]


def matrix_const_prod(c, a):
    return [[a[i][j] * c for i in range(len(a))] for j in range(len(a))]


def seidel(a, f, eps):
    print("Solving...")
    x = [0 for _ in f]
    n = len(x)
    converge = False
    iter_count = 0
    while not converge:
        x_new = [e for e in x]
        for i in range(n):
            s_new = sum(a[i][j] * x_new[j] for j in range(i))
            s = sum(a[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (f[i] - s_new - s) / a[i][i]

        converge = all(e < eps for e in map(lambda y, z: abs(y - z), x_new, x))
        iter_count += 1
        x = [e for e in x_new]
        if iter_count > 20000:
            break
    if converge:
        print(f"Successfully solved, iteration count: {iter_count}\n")
    else:
        print(f"!!! Unable to solve: n = {n}, alpha = {f[n - 1] - 1}")
    return x


def solve_fred2(x, n, k, phi, f, lmbd, eps):
    print("Calculating alpha-matrix...")
    a = alhpa_mat(phi, x, n)
    print("Done!")
    print("Calculating beta-matrix...")
    b = beta_mat(phi, k, x, n)
    print("Done!")
    print("Calculating gamma-vector...")
    g = gamma_vec(phi, k, f, x, n)
    print("Done!")
    mat = matrix_sub(a, matrix_const_prod(lmbd, b))
    print("Calculating coefficients c...")
    c = seidel(mat, g, eps)
    print("Done!")
    print("Calculating the solution...")
    u = [f(x[i]) + sum([c[j] * phi(j, x[i], x) for j in range(n)]) for i in range(n)]
    print("Done!")
    d = pd.DataFrame({"x": x, "u": u})
    d.to_csv(f"src/n{n}.csv")
    print("Dataframe saved!")
    if n > 8:
        print(trap_ints['8'])
        print(trap_ints['16'])
        print(f'Integral on {n} dots is better than the one on {int(n) / 2} by {trap_ints[str(int(int(n)/2))][0] - trap_ints[str((n))][0]}')
    return u


def seidel(a, f, eps):
    print("Solving...")
    x = [0 for _ in f]
    n = len(x)
    converge = False
    iter_count = 0
    while not converge:
        x_new = [e for e in x]
        for i in range(n):
            s_new = sum(a[i][j] * x_new[j] for j in range(i))
            s = sum(a[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (f[i] - s_new - s) / a[i][i]

        converge = all(e < eps for e in map(lambda y, z: abs(y - z), x_new, x))
        iter_count += 1
        x = [e for e in x_new]
        if iter_count > 20000:
            break
    if converge:
        print(f"Successfully solved, iteration count: {iter_count}\n")
    else:
        print(f"!!! Unable to solve: n = {n}, alpha = {f[n - 1] - 1}")
    return x


def min_sq(x_vals, y_vals, m):
    a = np.zeros((m + 1, m + 1))
    for i in range(m + 1):
        for j in range(m + 1):
            a[i][j] = sum([x ** (i + j) for x in x_vals])

    b = np.zeros(m + 1)
    for k in range(m + 1):
        b[k] = sum([x ** k * y for x, y in zip(x_vals, y_vals)])

    c = seidel(a, b, 0.01)
    return c


def poly_f(c, x):
    return sum([c[i] * x ** i for i in range(len(c))])




