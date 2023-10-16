from math import exp


def f(x):
    return exp(-x) * x


def g(x):
    return x - 2


def gol_sec_search(func, a, b, eps):
    iter_count = 1
    gol_rat = (5 ** 0.5 + 1) / 2  # golden-section ratio
    c = b - (b - a) / gol_rat
    d = a + (b - a) / gol_rat
    while abs(c - d) > eps:
        if func(c) < func(d):
            b = d
        else:
            a = c
        iter_count += 1
        c = b - (b - a) / gol_rat
        d = a + (b - a) / gol_rat
    return (c + d) / 2, iter_count


def penalty_function(func, cond, x, penalty):
    return func(x) + penalty * max(cond(x), 0) ** 2


def penalty_method(func, cond, penalty_func, init_penalty, a, b, eps, max_iter=1000, diff=5):
    penalties = [init_penalty]
    x_values = []
    min_x = None
    min_fx = float('inf')
    iter_s = 0
    for penalty in penalties:
        def pen_func(x):
            return penalty_func(func, cond, x, penalty)

        x_opt, iters = gol_sec_search(pen_func, a, b, eps)
        iter_s += iters
        fx = func(x_opt)
        if fx < min_fx:
            min_fx = fx
            min_x = x_opt

        if x_values:
            if abs(fx) < eps:
                break
        x_values.append(x_opt)

        if iter_s >= max_iter:
            print(f"Iteration count reached the limit of {max_iter}")
            return 0, max_iter

        penalties.append(penalty * diff)

    return min_x, iter_s



