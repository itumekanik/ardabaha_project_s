import numpy as np
np.set_printoptions(suppress=True, precision=2)

nodes = dict()
elements = dict()

class Node:
    def __init__(self, id, X, Y):
        self.id = id
        self.X = X
        self.Y = Y
        self.rest = [0, 0, 0]
        self.dof = [-1, -1, -1]
        self.load = [0, 0, 0]
        nodes[id] = self

class Frame2DElement:
    def __init__(self, id, n1, n2, E, b, h, mV=2500):
        self.id = id
        self.n1 = nodes[n1]
        self.n2 = nodes[n2]
        self.E = E
        self.b = b
        self.h = h
        self.A = b * h
        self.EA = E * b * h
        self.EIz = E * b * h**3 / 12
        self.mV = mV
        self.q = [0, 0, 0]  # [qn, qm, mz]
        elements[id] = self

    def code_vec(self):
        return [self.n1.dof[0], self.n1.dof[1], self.n1.dof[2],  
                self.n2.dof[0], self.n2.dof[1], self.n2.dof[2]]
    
    def L(self):
        X1 = self.n1.X
        X2 = self.n2.X
        Y1 = self.n1.Y
        Y2 = self.n2.Y
        return  ((X2-X1)**2 + (Y2-Y1)**2)**0.5

    def KL(self):
        EA = self.EA
        EIz = self.EIz
        L = self.L()
        return np.asarray([
            [EA/L,  0,            0,          -EA/L,  0,             0          ],
            [0,     12*EIz/L**3,  6*EIz/L**2,  0,    -12*EIz/L**3,   6*EIz/L**2 ],
            [0,     6*EIz/L**2,   4*EIz/L,     0,    -6*EIz/L**2,    2*EIz/L    ],
            [-EA/L, 0,            0,           EA/L,  0,             0          ],
            [0,    -12*EIz/L**3, -6*EIz/L**2,  0,     12*EIz/L**3,  -6*EIz/L**2 ],
            [0,     6*EIz/L**2,   2*EIz/L,     0,    -6*EIz/L**2,    4*EIz/L    ]
        ])
    
    def ML(self):
        mA = self.A * self.mV
        L = self.L()
        a = 0.5 * L
        return mA * a / 105 * np.asarray([
            [70,    0,            0,           35,     0,             0          ],
            [0,     78,           22*a,        0,      27,           -13*a       ],
            [0,     22*a,         8*a**2,      0,      13*a,         -6*a**2     ],
            [35,    0,            0,           70,     0,             0          ],
            [0,     27,           13*a,        0,      78,           -22*a       ],
            [0,    -13*a,        -6*a**2,      0,     -22*a,          8*a**2     ]
        ])


    def T(self):
        X1 = self.n1.X
        X2 = self.n2.X
        Y1 = self.n1.Y
        Y2 = self.n2.Y
        Lx = X2 - X1
        Ly = Y2 - Y1
        L = ((X2-X1)**2 + (Y2-Y1)**2)**0.5
        cos = Lx  / L
        sin = Ly / L
        return np.asarray([
            [cos, sin,  0,  0,   0,   0],
            [-sin, cos, 0,  0,   0,   0],
            [0,  0,     1,  0,   0,   0],
            [0,  0,     0,  cos, sin, 0],
            [0,  0,     0, -sin, cos, 0],
            [0,  0,     0,  0,   0,   1]
        ])

    def KG(self):
        return np.linalg.inv(self.T()) @ self.KL() @ self.T()
    
    def MG(self):
        return np.linalg.inv(self.T()) @ self.ML() @ self.T()

    def qL(self):
        X1 = self.n1.X
        X2 = self.n2.X
        Y1 = self.n1.Y
        Y2 = self.n2.Y
        L = ((X2-X1)**2 + (Y2-Y1)**2)**0.5
        qn, qm, mz = self.q  # unbox q vector in three different variables
        return np.asarray([
            qn*L/2, qm*L/2-mz, qm*L**2/12, qn*L/2, qm*L/2+mz, -qm*L**2/12
        ])

    def qG(self):
        return np.linalg.inv(self.T()) @ self.qL()



# Creating Node and Element objects
Node(id=1, X=0, Y=0)
Node(id=2, X=6, Y=0)
Node(id=3, X=0, Y=3)
Node(id=4, X=6, Y=3)

Frame2DElement(id=1, n1=1, n2=3, E=28e6, b=0.30, h=0.60)
Frame2DElement(id=2, n1=3, n2=4, E=28e6, b=0.30, h=0.60)
Frame2DElement(id=2, n1=4, n2=2, E=28e6, b=0.30, h=0.60)


# Defining restraint and load conditions 
nodes[1].rest = [1, 1, 1]
nodes[2].rest = [1, 1, 1]
nodes[3].load = [100, 0, 0]
elements[2].q = [0, -10, 0]  # Distributed load


# Dof Numbering
counter = 0
for i in range(3):
    for id, node in nodes.items():
        if node.rest[i] == 0:
            node.dof[i] = counter
            counter += 1

N = counter  # Number of unknown Dofs

for i in range(3):
    for id, node in nodes.items():
        if node.rest[i] == 1:
            node.dof[i] = counter
            counter += 1

M = counter # Number of All Dofs

print("Node Info")
for id, node in nodes.items():
    print(id, node.X, node.Y, node.rest, node.dof, node.load)

print("Unknown Dof Count (N):", N)
print("Total Dof Count   (M):", M)


# Stiffness Matrix Assembledge (KS)
KS = np.zeros((M, M))
for id, elm in elements.items():
    cv = elm.code_vec()
    KE = elm.KG()
    for i in range(len(cv)):
        for j in range(len(cv)):
            KS[cv[i], cv[j]] += KE[i, j]

# MASS Matrix Assembledge (MS)
MS = np.zeros((M, M))
for id, elm in elements.items():
    cv = elm.code_vec()
    ME = elm.ML()
    for i in range(len(cv)):
        for j in range(len(cv)):
            MS[cv[i], cv[j]] += ME[i, j]
            
    
# Distributed force vector Assembledge (qS)
qS = np.zeros(M)
for id, elm in elements.items():
    cv = elm.code_vec()
    qE = elm.qG()
    for i in range(len(cv)):
        qS[cv[i]] += qE[i]


# Concentrated force vector Assembledge (PS)
PS = np.zeros(M)
for id, node in nodes.items():
    cv = node.dof
    PE = node.load
    for i in range(len(cv)):
        PS[cv[i]] += PE[i]


# Equation solution
K11 = KS[0:N, 0:N]
K12 = KS[0:N, N:M]
K21 = KS[N:M, 0:N]
K22 = KS[N:M, N:M]

q1 = qS[0:N]
q2 = qS[N:M]

P1 = PS[0:N]


# u1 = np.linalg.inv(K11) @ (q1 + P1)
# P2 = K21 @ u1 - q2


# print("\n--- Displacements ---")
# for i, val in enumerate(u1):
#     print(f"u{i}={val}")

# print("\n--- Support Reactions ---")
# for i, val in enumerate(P2):
#     print(f"P{i+N}={val}")


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



f = interpolate.interp1d(all_data[:, 0], all_data[:, 1], fill_value=0, bounds_error=False)



dt = 0.001
t_array = np.arange(0, 60, dt)

# plt.plot(t_array, f(t_array), '-')
# plt.show()
# for t in t_array:
#     print(t)

M11 = MS[0:N, 0:N]
MI = np.linalg.inv(M11)

MIK = MI @ K11
C = 0.01 * MIK

unit_vec = np.zeros(N)
unit_vec[[0, 1]]=1.0
# unit_vec[:]=1.0

delta = [0]

u = np.zeros(N)
v = np.zeros(N)
for f in f(t_array):
    zdd = unit_vec * f
    a = -(zdd + C @ v + MIK @ u)
    v += dt * a
    u += dt * v
    delta.append(u[0])

# plt.plot(np.arange(0, 50, dt), delta[:], '-')
plt.plot(t_array, delta[:-1], '-')
plt.show()
    