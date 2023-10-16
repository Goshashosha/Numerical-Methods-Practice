import matplotlib.pyplot as plt
from external import *
from time import time


max_iter = 10000
sizes = (8, 10, 20, 50, 100, 500, 1000)
results = {'n': [], 'time': [], 'eigen_vec': []}

for n in sizes:
    start_time = time()
    eigen_val = rotations(gen_m(n), True)
    print(f"{n} : {eigen_val}")
    results['n'].append(n)
    results['time'].append(time() - start_time)
    results['eigen_vec'].append(eigen_val.tolist())

plt.plot(results['n'], results['time'])
plt.xlabel("Size of the matrix")
plt.ylabel("Calculation time")
plt.title("Dependence of calculation time on size of the matrix")
plt.savefig(f"src/all.png")
plt.show()


write_m(results['eigen_vec'], "src/eigen_vals")
