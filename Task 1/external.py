import sys


def norm_1(matrix):
    return max(sum(abs(matrix[i][j]) for i in range(len(matrix))) for j in range(len(matrix)))


def norm_inf(matrix):
    return max(sum(abs(matrix[i][j]) for j in range(len(matrix))) for i in range(len(matrix)))


def norm_euclid(matrix):
    sum(sum(matrix[i][j] ** 2 for j in range(len(matrix))) for i in range(len(matrix))) ** (1/2)


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

print(lu_decompose([[60, 91, 26], [60, 3, 75], [45, 90, 31]]))


