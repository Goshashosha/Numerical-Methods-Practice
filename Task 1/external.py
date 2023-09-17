import sys


def generate_matrix(n):
    print(f"Generating the {n}x{n} matrix")
    return [[min(i + 1, j + 1) for j in range(n)] for i in range(n)]


def matrix_multiplication(matrix1, matrix2):
    n = len(matrix1)
    res = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                res[i][j] += matrix1[i][k] * matrix2[k][j]
    return res


def lu_decompose(matrix):
    n = len(matrix)
    # implementation of gaussian algorithm
    u = [[matrix[i][j] for j in range(n)] for i in range(n)]
    l = [[int(i == j) for i in range(n)] for j in range(n)]
    for i in range(n):
        factors = [u[k + 1][i] / u[i][i] for k in range(n - 1)]
        if abs(u[i][i]) < sys.float_info.epsilon:
            print("The given matrix can not be LU-decomposed")
            exit()
        for j in range(i + 1, n):
            l[j][i] = factors[j - 1]
        for j in range(i + 1, n):
            u[j] = [u[j][k] - u[i][k] * factors[j - 1] for k in range(n)]
    return l, u


def invert(matrix):
    n = len(matrix)
    print("Decomposing...")
    l, u = lu_decompose(matrix)
    print("Decomposition complete!\n")
    print("Inverting...")
    l_inv = [[float(i == j) for i in range(n)] for j in range(n)]
    # inverting lower triangle matrix
    for i in range(n - 1):
        factors = [l[k][i] for k in range(i + 1, n)]
        for j in range(i + 1, n):
            l_inv[j] = [l_inv[j][k] - l_inv[i][k] * factors[j - i - 1] for k in range(n)]
            l[j] = [l[j][k] - l[i][k] * factors[j - i - 1] for k in range(n)]
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
    # calculation of inverted matrix
    m_inv = matrix_multiplication(u_inv, l_inv)
    print("Matrix inversion complete!")
    if len(matrix) <= 10:
        print("Inverse matrix:")
        for row in m_inv:
            print(' '.join(map(str, row)))
    return m_inv


def norm_1(matrix):
    return max(sum(abs(matrix[i][j]) for i in range(len(matrix))) for j in range(len(matrix)))


def norm_inf(matrix):
    return max(sum(abs(matrix[i][j]) for j in range(len(matrix))) for i in range(len(matrix)))


def norm_euclid(matrix):
    return sum(sum(matrix[i][j] ** 2 for j in range(len(matrix))) for i in range(len(matrix))) ** (1 / 2)


def condition_numbers(m):
    m_inv = invert(m)
    print("Calculating the condition numbers")
    return norm_1(m) * norm_1(m_inv), norm_inf(m) * norm_inf(m_inv), norm_euclid(m) * norm_euclid(m_inv)
