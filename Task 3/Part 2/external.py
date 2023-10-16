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
    main_vec = np.zeros(n)

    for i in range(n - 2):
        for k in range(i, n):
            main_vec[k] = a[n - 1, i] * a[n - 1, k]

        for j in range(n - 2, i, -1):
            main_vec[i] = main_vec[i] + a[j, i] * a[j, i]

        for k in range(i + 1, n):
            for j in range(n - 2, k - 1, -1):
                main_vec[k] = main_vec[k] + a[j, i] * a[j, k]
            for j in range(k - 1, i, -1):
                main_vec[k] = main_vec[k] + a[j, i] * a[k, j]

        coef_a = main_vec[i] ** (1/2)
        if abs(a[i + 1, i]) > eps:
            coef_b = 1 / coef_a

            for j in range(i + 2, n):
                a[j, i] = a[j, i] * coef_b

            main_vec[i], a[i + 1, i], g, main_vec2 = a[i + 1, i] * coef_b + np.sign(a[i + 1, i]), coef_a, 1 / abs(main_vec[i]), 0
            for k in range(i + 2, n):
                main_vec[k] = main_vec[k] * coef_b * g + np.sign(main_vec[i]) * abs(a[k, i + 1])
                main_vec2 = main_vec[k] * a[k, i] + main_vec2

            main_vec2 = g * main_vec2
            for k in range(i + 2, n):
                a[k, k] = a[k, k] - 2 * a[k, i] * main_vec[k] + main_vec2 * a[k, i] ** 2
                for j in range(k + 1, n):
                    a[j, k] = a[j, k] - a[j, i] * main_vec[k] - a[k, i] * main_vec[i] + main_vec2 * a[j, i] * a[k, i]
        else:
            if abs(coef_a) > eps:
                coef_b = 1 / coef_a
                for j in range(i + 2, n):
                    a[j, i] = a[j, i] * coef_b

                main_vec[i], a[i + 1, i], g, main_vec2 = -1, coef_a, 1, 0
                for k in range(i + 2, n):
                    main_vec[k] = main_vec[k] * coef_b * g + np.sign(main_vec[i]) * abs(a[k, i + 1])
                    main_vec2 = main_vec[k] * a[k, i] + main_vec2
                main_vec2 = g * main_vec2

                for k in range(i + 2, n):
                    a[k, k] = a[k, k] - 2 * a[k, i] * main_vec[k] + main_vec2 * a[k, i] ** 2
                    for j in range(k + 1, n):
                        a[j, k] = a[j, k] - a[j, i] * main_vec[k] - a[k, i] * main_vec[i] + main_vec2 * a[j, i] * a[k, i]
            else:
                main_vec[i] = 1

        return np.diag(a), np.diag(a, k=-1)
