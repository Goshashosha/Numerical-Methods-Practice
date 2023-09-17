from external import *
import matplotlib.pyplot as plt

iter_limit = 10000
n_limit = 100

for n in range(10, n_limit + 1, 10):
    a, f = generate_input(n)
    results = {'omega': [], 'iter_count': []}
    for i in range(21):
        omega = i / 10
        x, iter_count = alter_tri_method(a, f, omega, iter_limit)
        results['omega'].append(omega)
        results['iter_count'].append(iter_count)
    plt.plot(results['omega'], results['iter_count'], label=n)


plt.xlabel("Omega")
plt.ylabel("Iteration count")
plt.title("Dependence of iteration count on omega and n")
plt.legend()
plt.savefig(f"src/test_10-{n_limit}_{iter_limit}.png")
plt.show()
