from numpy import zeros, sin, pi, linspace, arcsin
import matplotlib.pyplot as plt


def stormer_method(f, y0, yp0, h, num_steps):
    y = zeros(num_steps)

    y[0] = y0
    y[1] = y0 + h * yp0

    for i in range(1, num_steps - 1):
        y[i + 1] = 2 * y[i] - y[i - 1] + (h ** 2 / 12) * (f(y[i + 1]) + 10 * f(y[i]) + f(y[i - 1]))

    return y


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
    return x


def min_sq(x_vals, y_vals, m):
    a = zeros((m + 1, m + 1))
    for i in range(m + 1):
        for j in range(m + 1):
            a[i][j] = sum([x ** (i + j) for x in x_vals])

    b = zeros(m + 1)
    for k in range(m + 1):
        b[k] = sum([x ** k * y for x, y in zip(x_vals, y_vals)])

    c = seidel(a, b, 0.001)
    return c


def poly_f(c, x):
    return sum([c[i] * x ** i for i in range(len(c))])


def f(y):
    return -sin(y)


# Начальные условия
y0 = 0
yp0 = 1
ys = []
hs = [1., 0.5, 0.1, 0.05, 0.01]

# Параметры метода
for h in hs:
    n = int((4 * pi) / h)
    print(f"При h = {h} количество точек в разбиении {n}")
    # Решение уравнения методом Штёрмера
    y_values = stormer_method(f, y0, yp0, h, n)
    ys.append(y_values[-1])
    # Визуализация результатов
    x_values = linspace(0, 4 * pi, n)
    x_interpolated = linspace(0, 4 * pi, 100)
    y_interpolated = poly_f(min_sq(x_values, y_values, 100), x_interpolated)
    plt.plot(x_values, y_values, label=n)
    # plt.plot(x_interpolated, y_interpolated, label="Интерполяция")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()

plt.show()
hs.reverse()
plt.plot([i * 0.1 for i in range(len(ys) - 1)], [abs(ys[i] - ys[i - 1]) for i in range(1, len(ys))], label="Convergence")
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()