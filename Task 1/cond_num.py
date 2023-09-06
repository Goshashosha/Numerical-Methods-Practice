from time import time
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from external import condition_numbers, generate_matrix
from pandas import DataFrame


results = {'step': [], '1': [], 'inf': [], 'e': []}
timings = []
minute_limit = 1

for n in range(10, 10000, 10):
    start_time = time()
    m = generate_matrix(n)
    c_1, c_inf, c_e = condition_numbers(m)
    results['step'].append(n)
    results['1'].append(c_1)
    results['inf'].append(c_inf)
    results['e'].append(c_e)
    timings.append(time() - start_time)
    print(f"Condition number, norm 1 = {c_1}\n"
          f"Condition number, norm inf = {c_inf}\n"
          f"Condition number, norm e = {c_e}\n"
          f"Time spent = " + "{:.2f}".format(timings[-1]) + '\n')
    if timings[-1] >= 60 * minute_limit:
        break

df = DataFrame(results)
df.to_csv("results.csv", index=False)

x = [10 * (i  +  1) for i in range(len(timings))]
plt.plot(x, timings)
print(timings)
plt.show()