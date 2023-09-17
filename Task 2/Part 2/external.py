def generate_input(n):
    a = [[1 / (i + j + 1) for j in range(n)] for i in range(n)]
    return a, [sum(a[i]) for i in range(n)]


def mat_mat_mult(matrix1, matrix2):
    n = len(matrix1)
    res = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                res[i][j] += matrix1[i][k] * matrix2[k][j]
    return res


def lower_tri_invert(l):
    n = len(l)
    l_inv = [[float(i == j) for i in range(n)] for j in range(n)]
    # inverting lower triangle matrix
    for i in range(n):
        divider = l[i][i]
        l[i] = [l[i][j] / divider for j in range(n)]
        l_inv[i] = [l_inv[i][j] / divider for j in range(n)]
        factors = [l[k][i] for k in range(i + 1, n)]
        for j in range(i + 1, n):
            l_inv[j] = [l_inv[j][k] - l_inv[i][k] * factors[j - i - 1] for k in range(n)]
            l[j] = [l[j][k] - l[i][k] * factors[j - i - 1] for k in range(n)]
    return l_inv


def upper_tri_invert(u):
    n = len(u)
    u_inv = [[float(i == j) for i in range(n)] for j in range(n)]
    # inverting upper triangle matrix (reverse, because of its structure)
    for i in range(n - 1, -1, -1):
        # for the diagonal elem to be equal to 1
        divider = u[i][i]
        u[i] = [u[i][j] / divider for j in range(n)]
        u_inv[i] = [u_inv[i][j] / divider for j in range(n)]
        factors = [u[k][i] for k in range(i - 1, -1, -1)]
        for j in range(i - 1, -1, -1):
            u_inv[j] = [u_inv[j][k] - u_inv[i][k] * factors[j - i + 1] for k in range(n)]
            u[j] = [u[j][k] - u[i][k] * factors[j - i + 1] for k in range(n)]
    return u_inv


def mat_vec_mult(m, v):
    return [sum(m[i][j] * v[j] for j in range(len(v))) for i in range(len(v))]


def vec_subtr(v1, v2):
    return [v1[i] - v2[i] for i in range(len(v1))]


def vec_add(v1, v2):
    return [v1[i] + v2[i] for i in range(len(v1))]


def scal_prod(v1, v2):
    return sum(v1[i] * v2[i] for i in range(len(v1)))


def const_prod(c, v1):
    return [c * v1[i] for i in range(len(v1))]


def alter_tri_method(a, f, omega, iter_limit):
    n = len(f)
    y = [0 for _ in f]

    r = [[0 for _ in range(n)] for _ in range(n)]
    r_t = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            r[i][j] = a[i][j]
            if i == j:
                r[i][j] /= 2
        for j in range(i, n):
            r_t[i][j] = a[i][j]
            if i == j:
                r_t[i][j] /= 2

    # B = gamma_1 * gamma_2
    gamma_1 = [[omega * r_t[i][j] + (i == j) for j in range(n)] for i in range(n)]
    gamma_2 = [[omega * r[i][j] + (i == j) for j in range(n)] for i in range(n)]

    b_inv = mat_mat_mult(lower_tri_invert(gamma_2), upper_tri_invert(gamma_1))

    converge = False
    iter_count = 0
    while not converge:
        rk = vec_subtr(mat_vec_mult(a, y), f)
        wk = mat_vec_mult(b_inv, rk)
        tau = scal_prod(mat_vec_mult(a, wk), wk) / \
            scal_prod(mat_vec_mult(mat_mat_mult(b_inv, a), wk), mat_vec_mult(a, wk))

        y_new = vec_add(vec_subtr(y, const_prod(tau, mat_vec_mult(mat_mat_mult(b_inv, a), y))),
                        const_prod(tau, mat_vec_mult(b_inv, f)))

        iter_count += 1
        converge = all(e < 1e-1 for e in map(lambda x, z: abs(x - z), y_new, y))
        y = [e for e in y_new]
        if iter_count > iter_limit:
            break
    return y, iter_count
