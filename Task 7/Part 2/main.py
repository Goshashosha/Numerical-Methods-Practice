from matplotlib import pyplot as plt
from external import quad_solve_f1, linspace, k, f, poly_f, min_sq, reg_tikh, auto_k, auto_f, spline
from inteq import SolveFredholm


for eps in (1e-3, 1e-5, 1e-10):
    x, h = linspace(-1, 1, 30)
    u = quad_solve_f1(x, k, f, eps)
    # auto_u = SolveFredholm(auto_k, auto_f, gamma=eps, num=30)
    u_r = reg_tikh(-1, 1, 30, f, k, eps)
    x_spline, _ = linspace(-1, 1, 70)
    u_spline = spline(x, h, u, x_spline)
    u_r_spline = spline(x, h, u_r, x_spline)

    plt.plot(x, u_r, label="reg")
    # plt.plot(auto_u[0], auto_u[1], label="auto")
    plt.plot(x, u, label="nat")
    plt.plot(x_spline, u_spline, label="nat spline", linestyle='dashed')
    plt.plot(x_spline, u_r_spline, label="reg spline", linestyle='dashed')
    plt.xlabel("X")
    plt.ylabel("U")
    plt.title("Estimated solution U(x)")
    plt.legend()
    plt.show()
    plt.savefig(f"src/fig{int(1 / eps)}.png")
