from external import f, g, penalty_method, penalty_function
import matplotlib.pyplot as plt

ab = [float(e) for e in input("Please, enter the a and b values (through space) ").split()]
res, it_cnt = penalty_method(f, g, penalty_function, 10, ab[0], ab[1], 0.001, 100000, 10)
print(f"The minimum of f, as x >= 2, on [{ab[0]}, {ab[1]}] is found in point "
      f"{res}\nIteration count = {it_cnt}\n")

epss = [1 / 10 ** i for i in range(6, 15)]
it_cnts = []
for eps in epss:
    _, it_cnt = penalty_method(f, g, penalty_function, 10, ab[0], ab[1], eps, 100000, 10)
    it_cnts.append(it_cnt)
    print(f"eps = {eps}, it_cnt = {it_cnt}")

plt.plot(epss, it_cnts)
plt.xlabel("Tolerance")
plt.ylabel("Iteration count")
plt.title("Dependence of iteration count on tolerance value")
plt.savefig(f"src/Dependence.png")
plt.show()
