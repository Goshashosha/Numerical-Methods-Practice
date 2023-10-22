from math import cos, acos, pi


def cheb_f(x, n):
    return 2 ** (1 - n) * cos(n * acos(x))


def cheb_zeros(n):
    return [cos((2 * k + 1) * pi / 2 / n) for k in range(n)]


def runge_f(x):
    return 1 / (1 + 25 * x ** 2)


def linsps(a, b, n):
    return [a + (b - a) / n * k for k in range(0, n + 1)]


def save_poly(polynoms, cheb, n, x):
    name = "cheb" if cheb else "uni"
    with open(f"src/{name}-{n}.txt", "a") as f:
        f.write(f"For x = {x}:\n")
        for e in polynoms:
            f.write(f"i = {e[0]}, coefficient = {e[1]}\n")
        f.write("\n")


def lagr_basis_polynom(x_vals, i, x):
    res = 1
    for j in range(len(x_vals)):
        if j != i:
            res *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
    return res


def lagrange_interpolation(x_vals, y_vals, x, cheb=False):
    res = 0
    polynoms = []
    for i in range(len(x_vals)):
        lbp = lagr_basis_polynom(x_vals, i, x)
        polynoms.append([lbp, i])
        res += y_vals[i] * lbp

    save_poly(polynoms, cheb, len(x_vals), x)
    return res


if __name__ == "__main__":
    print(lagrange_interpolation([-1, 0, 1], [1, 0, 1], 2))
