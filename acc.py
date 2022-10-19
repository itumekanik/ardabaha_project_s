from _timing import timeit
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import acc_datas

import math
import time
import concurrent.futures

start = time.time()

data_path = r"acc_datas"
data_paths = [os.path.join(data_path, name) for name in os.listdir(data_path)]
all_data = []

for file in data_paths:
    data = []
    with open(file, "r") as f:
        all_data = np.asarray(
            [(lambda line: list(map(float, line.split("\t"))))(line)
             for line in f.readlines()],
            dtype=float)

all_data[:, 1] *= 9.81

dt = 0.0001
t_max = 30

f = interpolate.interp1d(all_data[:, 0], all_data[:, 1])
t_array = np.arange(0, t_max, dt)
zdd_array = f(t_array)



def S_acc (T: float):

    timeit["A1"]

    global max_acc

    v0 = 0.0
    x0 = 0.0
    w = 2.0 * np.pi / T
    w2 = w ** 2

    v_all = [v0]
    x_all = [x0]

    i = -1
    t0 = 0
    t_all = [t0]

    timeit["A1"]

    while t0 < 40:

        timeit["A2"]
        i += 1
        t0 = i * dt
        t1 = t0 + dt

        if t0 < t_max:
            zdd0 = f(t0)
            zdd1 = f(t1)
        else:
            zdd0 = 0
            zdd1 = 0

        timeit["A2"]


        timeit["A3"]

        zdd_ort = 0.5 * (zdd0 + zdd1)
        v1 = v0 - dt * (zdd_ort + w2 * x0 + 0.1 * w * v0)
        v_ort = 0.5 * (v0 + v1)
        x1 = x0 + dt * v_ort

        timeit["A3"]


        timeit["A4"]

        v_all.append(v1)
        x_all.append(x1)
        t_all.append(t1)

        timeit["A4"]

        v0 = v1
        x0 = x1

    # plt.plot(t_all, x_all, '-')
    # plt.show()

    max_acc = abs(max(x_all, key=abs))

    return max_acc

minT = 1
maxT = 5
spec_interval = 1

spec_acc = []
spec_t = (np.arange(minT,maxT,spec_interval))

# def do_something(seconds):
#     print(f'Sleeping {seconds} second(s)...')
#     time.sleep(seconds)
#     return 'Done Sleeping...'
#
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     secs = [5,4,3,2,1]
#     results = executor.map(do_something, secs)
#
#     for result in results:
#         print(result)


# with concurrent.futures.ProcessPoolExecutor() as executor:
#     results = executor.map(S_acc, spec_t)
#
#     for result in results:
#         print(result)
#
# print(spec_acc)


i=2
for t in np.arange(minT,maxT,spec_interval):
    perc = (i/len(spec_t))*100
    i += 1
    print("Process: %", perc)
    S_acc(t)
    spec_acc.append(max_acc)


d = interpolate.interp1d(spec_t, spec_acc)
spec_t_array = np.arange(minT, maxT-1, dt)
spec_array = d(spec_t_array)

end = time.time()
print("Toplam süre:", int((end - start)//60), "dk", int((end-start)%60), "sn" )

timeit.report()

plt.plot(spec_t_array, spec_array, '-')
plt.show()

