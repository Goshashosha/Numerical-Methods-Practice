from math import cos


def y_f(x):
    return 1 - cos(x)


def linspace(n):
    return [i / n for i in range(n)]


def zero_v(n):
    return [0. for _ in range(n)]


def zero_m(n):
    return [zero_v(n) for _ in range(n)]


def save_c(c, m):
    with open(f"src/coefficients-{m}.txt", "w") as f:
        for e in c:
            f.write(f"{e} ")


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
    a = zero_m(m + 1)
    for i in range(m + 1):
        for j in range(m + 1):
            a[i][j] = sum([x ** (i + j) for x in x_vals])

    b = zero_v(m + 1)
    for k in range(m + 1):
        b[k] = sum([x ** k * y for x, y in zip(x_vals, y_vals)])

    c = seidel(a, b, 0.01)
    save_c(c, m)
    return c


def poly_f(c, x):
    return sum([c[i] * x ** i for i in range(len(c))])
