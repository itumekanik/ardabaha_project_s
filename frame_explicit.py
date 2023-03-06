import numpy as np
np.set_printoptions(suppress=True, precision=0)

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
nodes[3].rest = [1, 1, 1]
nodes[2].rest = [1, 0, 0]
nodes[2].load = [0, -70, 0]
elements[1].q = [0, -15, 0]  # Distributed load
elements[2].q = [0, 0, -10]  # Distributed load


# Dof Numbering
counter = 0
for id, node in nodes.items():
    if node.rest[0] == 0:
        node.dof[0] = counter
        counter += 1
    if node.rest[1] == 0:
        node.dof[1] = counter
        counter += 1
    if node.rest[2] == 0:
        node.dof[2] = counter
        counter += 1

N = counter  # Number of unknown Dofs

for id, node in nodes.items():
    if node.rest[0] == 1:
        node.dof[0] = counter
        counter += 1
    if node.rest[1] == 1:
        node.dof[1] = counter
        counter += 1
    if node.rest[2] == 1:
        node.dof[2] = counter
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


u1 = np.linalg.inv(K11) @ (q1 + P1)
P2 = K21 @ u1 - q2


print("\n--- Displacements ---")
for i, val in enumerate(u1):
    print(f"u{i}={val}")

print("\n--- Support Reactions ---")
for i, val in enumerate(P2):
    print(f"P{i+N}={val}")


print(K11)
print(q1)
print(P1)