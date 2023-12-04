from matplotlib import pyplot as plt
from external import solve_fred2, linspace, k, f, phi, lbd, min_sq, poly_f, trap_int_2


for n in (8, 16, 32, 64):
    x = linspace(0, 1, n)
    u = solve_fred2(x, n, k, phi, f, lbd, 1e-4)
    c = min_sq(x, u, n)

    # checking the solution
    print([u[i] - trap_int_2(x, lambda s: abs(x[i] - s) * u[int(s * n)]) - f(x[i]) for i in range(n - 1)])

    plt.plot(x, u, label="nat")
    plt.plot(x, [poly_f(c, xx) for xx in x], label="spline")
    plt.xlabel("X")
    plt.ylabel("U")
    plt.title("Estimated solution U(x)")
    plt.legend()
    plt.show()
    plt.savefig(f"src/fig{n}.png")


