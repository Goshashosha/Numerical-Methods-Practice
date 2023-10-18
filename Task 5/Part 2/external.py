def f(x):
    return 10 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2


def norm(x):
    return sum(e ** 2 for e in (x))


def gol_sec_search(func, a, b, x, grad, eps):
    iter_count = 1
    gol_rat = (5 ** 0.5 + 1) / 2  # golden-section ratio
    c = b - (b - a) / gol_rat
    d = a + (b - a) / gol_rat
    while abs(c - d) > eps:
        if func([x[i] - c * grad[i] for i in range(2)]) < func([x[i] - d * grad[i] for i in range(2)]):
            b = d
        else:
            a = c
        iter_count += 1
        c = b - (b - a) / gol_rat
        d = a + (b - a) / gol_rat
    return (c + d) / 2, iter_count


def grad_f(x):
    df_dx0 = 40 * x[0]**3 - 40 * x[0] * x[1] + 2 * x[0] - 2
    df_dx1 = 20 * x[1] - 20 * x[0]**2
    return [df_dx0, df_dx1]


def grad_desc(func, grad_func, x0, eps=1e-3, max_iter=1e9):
    x = x0
    iter_s = 0

    while True:
        grad = grad_func(x)
        alpha, iter_cnt = gol_sec_search(func, 0, 1, x, grad, eps)

        print(x, alpha, grad)
        x_new = [x[i] - alpha * grad[i] for i in range(2)]

        if all([abs(x_new[i] - x[i]) < eps for i in range(2)]):
            break

        x = x_new
        iter_s += iter_cnt

        if iter_s >= max_iter:
            print(f"Iteration count reached the limit of {max_iter}")
            return 0, max_iter

    return x_new, iter_s


