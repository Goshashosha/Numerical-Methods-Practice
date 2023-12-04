import math
import numpy
import numpy as np
import pandas as pd


def auto_f(x):
    print(x)
    res = [2 * xx / 5 * (math.atan(xx / 5) - math.atan(1 / 5)) + math.log(26 / (xx ** 2 + 25), math.e) for xx in x]
    print(res)
    return res


def auto_k(x, s):
    return np.abs(x - s)


def k(x, s):
    return abs(x - s)


def f(x):
    return 2 * x / 5 * (math.atan(x / 5) - math.atan(1 / 5)) + math.log(26 / (x ** 2 + 25), math.e)


def linspace(a, b, n):
    return [a + i * (b - a) / (n - 1) for i in range(n)], (b - a) / (n - 1)


def con_grad_solve(x0, a, b, eps, max_iter=10000):
    x = x0
    r = b - np.dot(a, x)
    p = r
    rs_old = np.dot(np.transpose(r), r)

    for i in range(max_iter):
        ap = np.dot(a, p)
        alpha = rs_old / np.dot(np.transpose(p), ap)
        x = x + alpha * p
        r = r - alpha * ap
        rs_new = np.dot(np.transpose(r), r)

        if np.sqrt(rs_new) < eps:
            break

        p = r + (rs_new / rs_old) * p
        rs_old = rs_new
    print(x)

    return x


def quad_solve_f1(xs, k, f, eps):
    h = xs[1] - xs[0]
    a = [[h * k(x, s) for s in xs] for x in xs]
    b = [f(x) for x in xs]

    y = con_grad_solve([0 for _ in xs], a, b, eps)
    df = pd.DataFrame({'x': xs, 'y': y})
    df.to_csv(f"src/eps{int(1 / eps)}.csv")

    return y


def norm(x):
    return sum(e ** 2 for e in x)


def gol_sec_search(func, a, b, x, grad, eps):
    iter_count = 1
    gol_rat = (5 ** 0.5 + 1) / 2  # golden-section ratio
    c = b - (b - a) / gol_rat
    d = a + (b - a) / gol_rat
    while abs(c - d) > eps:
        if func([x[i] - c * grad[i] for i in range(x)]) < func([x[i] - d * grad[i] for i in range(x)]):
            b = d
        else:
            a = c
        iter_count += 1
        c = b - (b - a) / gol_rat
        d = a + (b - a) / gol_rat
    return (c + d) / 2, iter_count


def seidel(a, f, eps=1e-3):
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


def int_to_mat(a, b, n, k, f):
    x, h = linspace(a, b, n)
    s, _ = linspace(a, b, n)
    mat = numpy.zeros((n, n))
    vec = numpy.zeros(n)

    for i in range(n):
        for j in range(n):
            mat[i, j] = k(x[i], s[j]) * h
        vec[i] = f(x[i])

    return mat, vec


def find_alpha():
    return 0.0001


def reg_tikh(a, b, n, f, k, eps):
    mat, vec = int_to_mat(a, b, n, k, f)
    mat = np.array(mat)
    vec = np.array(vec)
    alpha = find_alpha()
    alpha_mat = np.diag(np.full(n, alpha))
    I = np.eye(n)
    return seidel(mat.T @ mat + alpha_mat @ I, mat.T @ vec, eps)


def get_mat_vec(x, h, y):
    n = len(y)
    a = np.zeros((n, n))
    b = np.zeros(n)

    # edge conditions:
    a[0, 0] = 2
    a[0, 1] = 1
    a[n - 1, n - 2] = 1
    a[n - 1, n - 1] = 2
    b[0] = 3 * (y[1] - y[0]) / (x[1] - x[0])
    b[n - 1] = 3 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])

    # constants (разбиение равномерно)
    lbd = 1 / 2
    mu = 1 / 2
    for i in range(1, n - 1):
        b[i] = 3 * (mu * (y[i + 1] - y[i]) / (x[i + 1] - x[i]) + lbd * (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
        a[i, i - 1] = lbd
        a[i, i] = 2
        a[i, i + 1] = mu

    return a, b


def solve_tridiagonal(a, b):
    n = len(b)
    # Прямой ход
    for i in range(1, n):
        m = a[i][0] / a[i - 1][1]
        a[i][1] -= m * a[i - 1][2]
        b[i] -= m * b[i - 1]

    # Обратный ход
    x = np.zeros(n)
    x[-1] = b[-1] / a[-1][1]
    for i in range(n - 2, -1, -1):
        x[i] = (b[i] - a[i][2] * x[i + 1]) / a[i][1]

    return x


def interpolate(x, y, m, x_spline):
    res = []
    for x_s in x_spline:
        i = np.searchsorted(x, x_s) - 1
        if i < 0:
            i = 0
        print(x_spline)
        print(x)
        print(y)
        print(m)
        print(i)
        res.append(y[i] / (x[i + 1] - x[i]) ** 3 * (x[i + 1] - x_s) ** 2 * (x[i + 1] + 2 * x_s - 3 * x[i]) +
                   y[i + 1] / (x[i + 1] - x[i]) ** 3 * (x_s - x[i]) ** 2 * (3 * x[i + 1] - 2 * x_s - x[i]) +
                   m[i] / (x[i + 1] - x[i]) ** 2 * (x[i + 1] - x_s) ** 2 * (x_s - x[i]) +
                   m[i + 1] / (x[i + 1] - x[i]) ** 2 * (x_s - x[i]) ** 2 * (x_s - x[i + 1]))
    return res


def spline(x, h, y, x_spline):
    a, b = get_mat_vec(x, h, y)
    m = seidel(a, b, eps=1e-10)
    return interpolate(x, y, m, x_spline)
