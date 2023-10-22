from external import *
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(13, 7))
n = 10

for m in (1, 2, 3, 4):
    x = linspace(n)
    y = [y_f(x[i]) for i in range(n)]

    coefficients = min_sq(x, y, m)
    y_poly = [poly_f(coefficients, x[i]) for i in range(n)]

    plt.subplot(1, 4, m)
    plt.plot(x, y_poly, 'r')
    plt.plot(x, y, linestyle="dashed", color='g')
    plt.title(f"m = {m}")

plt.show()
# plt.savefig(f"src/plot.png")
