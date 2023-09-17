import matplotlib.pyplot as plt
from external import *
from pandas import DataFrame


iter_limit = 20000
eps = 1e-3
converge_issue = {'n': [], 'alpha': []}

for n in range(10, 41, 10):
    results = {'alpha': [], 'iter_count': []}
    alpha = 0
    for _ in range(100):
        alpha += 0.01
        print(alpha, "alpha")
        a, f = generate_input(n, alpha)
        x, iter_count = seidel(a, f, eps)
        if any(abs(e) > eps for e in map(lambda y: y - 1, x)):
            converge_issue['n'].append(n)
            converge_issue['alpha'].append(alpha)
        results['alpha'].append(alpha)
        results['iter_count'].append(iter_count)
    plt.plot(results['alpha'], results['iter_count'], label=n)

df = DataFrame(converge_issue)
df.to_csv(f"src/error_{iter_limit}.csv", index=False)

plt.xlabel("Alpha")
plt.ylabel("Iteration count")
plt.title("Dependence of iteration count on alpha and n")
plt.legend()
plt.savefig(f"src/test_10-40_{iter_limit}.png")
plt.show()
