import numpy as np
import matplotlib.pyplot as plt


# constants
l = 1
mu1 = 0
mu2 = 0


# the given f (right side of the equation)
def f(t, u):
    return -np.exp(u)


# Seidel's method for solving SLAE
def seidel(a, b, eps):
    print("Solving...")
    x = [0 for _ in b]
    n = len(x)
    converge = False
    iter_count = 0
    while not converge:
        x_new = [e for e in x]
        for i in range(n):
            s_new = sum(a[i][j] * x_new[j] for j in range(i))
            s = sum(a[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s_new - s) / a[i][i]

        converge = all(e < eps for e in map(lambda y, z: abs(y - z), x_new, x))
        iter_count += 1
        x = [e for e in x_new]
        if iter_count > 20000:
            break
    if converge:
        print(f"Successfully solved, iteration count: {iter_count}\n")
    else:
        print(f"!!! Unable to solve: n = {n}, alpha = {b[n - 1] - 1}")

    d = x[0]
    x = [xx - d for xx in x]
    return x


# Norm in R^2
def norm(xs):
    return sum(x ** 2 for x in xs) ** (1 / 2)


# Approximate derivative of f
def f_der(x, u_prev, u, h):
    return (f(x, u) - f(x, u_prev)) / h


# Implementation of Newton's method
def newton_method(m1, m2, lower_bound, upper_bound, n, eps):
    x_vals, h = np.linspace(lower_bound, upper_bound, n, retstep=True)
    u = np.zeros(n)
    u[0] = m1
    u[-1] = m2
    u_prev = np.zeros(n)
    u_prev[1] = 1
    ic = 0
    while norm(u[1:-1] - u_prev[1:-1]) > eps:
        u[0] = m1
        u[-1] = m2
        u_prev = np.copy(u)
        mat = np.zeros((n, n))
        vec = np.zeros(n)

        for i in range(n):
            for j in range(n):
                if i == j:
                    if i == n - 1:
                        mat[i, j] = -2 / (h ** 2) - f_der(x_vals[i], u_prev[i - 1], u_prev[i], h)
                    else:
                        mat[i, j] = -2 / (h ** 2) - f_der(x_vals[i], u_prev[i], u_prev[i + 1], h)
                elif abs(i - j) == 1:
                    mat[i, j] = 1 / (h ** 2)

        for i in range(n):
            if i == n - 1:
                vec[i] = f(x_vals[i], u_prev[i]) - f_der(x_vals[i], u_prev[i - 1], u_prev[i], h) * u_prev[i]
            else:
                vec[i] = f(x_vals[i], u_prev[i]) - f_der(x_vals[i], u_prev[i], u_prev[i + 1], h) * u_prev[i]

        u = seidel(mat, vec, eps)

        print(u_prev, u)
        # plt.plot(x_vals, u, label=ic)

        ic += 1
        if ic > 10000:
            break

    print(ic)
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.legend()
    # plt.savefig(f"src/figx{n}.png")
    # plt.show()
    return x_vals, u, ic


res = []
epss = (1e-2, 1e-3, 1e-4, 1e-5)
for eps in epss:
    x, y, it_cnt = newton_method(mu1, mu2, 0, 1, 50, eps)
    res.append(it_cnt)

plt.plot(epss, res)
plt.show()
plt.savefig(f"src/convergence-rate.png")
