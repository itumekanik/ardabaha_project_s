import numpy as np

zdd = np.asarray([3,4,5,6,7,8], dtype=float)



print(np.roll(zdd, -1) + zdd)