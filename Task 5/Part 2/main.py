from external import grad_desc, f, grad_f
import matplotlib.pyplot as plt


x0 = [float(e) for e in input("Please, enter the x0 value (through space) ").split()]
res, it_cnt = grad_desc(f, grad_f, x0)
print(f"The minimum of f is found in point "
      f"{res}\nIteration count = {it_cnt}\n")

exit()

epss = [1 / 10 ** i for i in range(6, 15)]
it_cnts = []
for eps in epss:
    _, it_cnt = grad_desc(f, grad_f, [0.8, 0.6], eps, 1e9)
    it_cnts.append(it_cnt)
    print(f"eps = {eps}, it_cnt = {it_cnt}")

plt.plot(epss, it_cnts)
plt.xlabel("Tolerance")
plt.ylabel("Iteration count")
plt.title("Dependence of iteration count on tolerance value")
plt.savefig(f"src/Dependence.png")
plt.show()