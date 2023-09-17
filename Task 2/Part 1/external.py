def generate_input(n, alpha):
    print(f"Generating A and f, size: {n}...")
    a = [[0 for _ in range(n)] for _ in range(n)]
    a[0][0], a[0][1] = 2, -1 - alpha
    a[n - 1][n - 2], a[n - 1][n - 1] = -1 + alpha, 2
    for i in range(1, n - 1):
        for j in range(i - 1, i + 2):
            if i == j:
                a[i][j] = 2
            elif i < j:
                a[i][j] = -1 - alpha
            else:
                a[i][j] = -1 + alpha
    f = [0 for _ in range(n)]
    f[0], f[n - 1] = 1 - alpha, 1 + alpha
    print("Generating process complete!\n")
    return a, f


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
    return x, iter_count