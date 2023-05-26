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



f = interpolate.interp1d(all_data[:, 0], all_data[:, 1])
# t_array = np.arange(0, t_max, dt)
# zdd_array = f(t_array)



def S_acc (T: float):

    # timeit["A1"]

    global max_acc

    v0 = 0.0
    x0 = 0.0
    w = 2.0 * np.pi / T
    w2 = w ** 2

    a_all = [0]
    v_all = [v0]
    x_all = [x0]

    i = -1
    t0 = 0
    t_all = [t0]

    # timeit["A1"]

    dt = 0.0001 * 10 * T
    t_max = 31      

    t = np.arange(0, t_max, dt)

    zdd = f(t)

    zdd_ort_array =  0.5 * (np.roll(zdd, -1) + zdd)

    zdd_max = abs(max(zdd_ort_array, key=abs)) # Maksimum yer ivmesi



    for i, t0 in enumerate(t[:-1]):
        t1 = t[i+1]

        # timeit["A3"]

        # zdd_ort = 0.5 * (zdd0 + zdd1)
        zdd_ort = zdd_ort_array[i]
        v1 = v0 - dt * (zdd_ort + w2 * x0 + 0.05 * w * v0)
        v_ort = 0.5 * (v0 + v1)
        x1 = x0 + dt * v_ort

        # timeit["A3"]


        # timeit["A4"]

        a_all.append(zdd_ort + w2 * x0 + 0.05 * w * v0)
        v_all.append(v1)
        x_all.append(x1)
        t_all.append(t1)

        # timeit["A4"]

        v0 = v1
        x0 = x1

    # plt.plot(t_all, x_all, '-')
    # plt.show()

    max_acc = abs(max(a_all, key=abs))/zdd_max

    return max_acc

minT = 0.01
maxT = 8
spec_interval = 0.02

spec_acc = []
spec_t = (np.arange(minT,maxT,spec_interval))


i=1
for t in np.arange(minT,maxT,spec_interval):
    perc = (i/len(spec_t))*100
    i += 1
    print("Process: %", perc)
    S_acc(t)
    spec_acc.append(max_acc)


d = interpolate.interp1d(spec_t, spec_acc)
spec_t_array = np.arange(minT, maxT-1, 0.0001)
spec_array = d(spec_t_array)

end = time.time()
print("Toplam süre:", int((end - start)//60), "dk", (end-start)%60, "sn" )

# timeit.report()

plt.plot(spec_t_array, spec_array, '-')
plt.show()

