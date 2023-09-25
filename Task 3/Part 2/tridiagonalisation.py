import matplotlib.pyplot as plt
from external import *
from time import time


results = {'n': [], 'time': []}

for n in range(10, 1001, 10):
    start_time = time()
    main_d, sec_d = hh_trans(gen_m(n))
    print(f"{n} : {main_d}\n{sec_d}\n")
    results['n'].append(n)
    results['time'].append(time() - start_time)

plt.plot(results['n'], results['time'])
plt.xlabel("Size of the matrix")
plt.ylabel("Calculation time")
plt.title("Dependence of calculation time on size of the matrix")
plt.savefig(f"src/all.png")
plt.show()
