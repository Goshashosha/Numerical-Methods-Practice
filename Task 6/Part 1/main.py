from external import *
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(12, 10))

for n, i in [(4, 1), (6, 3), (10, 5), (20, 7)]:
    x_uni = linsps(-1, 1, n)
    x_cheb = cheb_zeros(n)
    y_uni = [runge_f(x) for x in x_uni]
    y_cheb = [runge_f(x) for x in x_cheb]

    y_int_uni = [lagrange_interpolation(x_uni, y_uni, x) for x in x_uni]
    y_int_cheb = [lagrange_interpolation(x_cheb, y_cheb, x, True) for x in x_cheb]
    y_runge = [runge_f(x) for x in x_uni]

    lsp = [-1 + i / 100 for i in range(201)]

    plt.subplot(4, 2, i)
    plt.plot(x_uni, y_int_uni, 'r')
    plt.plot(lsp, [runge_f(x) for x in lsp], linestyle="dashed", color='g')
    if i == 1:
        plt.title("Uniform")

    plt.subplot(4, 2, i + 1)
    plt.plot(x_cheb, y_int_cheb, 'r')
    plt.plot(lsp, [runge_f(x) for x in lsp], linestyle="dashed", color='g')
    if i == 1:
        plt.title("Chebyshov")

plt.legend()
plt.show()
plt.savefig(f"src/plot.png")
