#import scipy as sc
import numpy as np
import matplotlib.pyplot as plt

def FDLaplacian(NX, NY):
    Dx = np.zeros((NX, NX - 1))
    Dy = np.zeros((NY, NY - 1))
    np.fill_diagonal(Dx, 1)
    np.fill_diagonal(Dy, 1)
    np.fill_diagonal(Dx[1:, :], -1)
    np.fill_diagonal(Dy[1:, :], -1)
    Lx = np.dot(Dx.T, Dx)
    Ly = np.dot(Dy.T, Dy)
    Ax = np.kron( np.eye(np.shape(Ly)[1]), Lx)
    Ay = np.kron(Ly, np.eye(np.shape(Lx)[1]))
    return Ax + Ay

def f(x, y):
    value = 0
    alpha = 40
    for i in range(1, 10):
        for j in range(1, 5):
            value += np.exp(-alpha *(x - i) ** 2 - alpha * (y - j) ** 2)
    return value

#globals
NX = 4
NY = 4
xdim = 10
ydim = 5
h1 = xdim / NX
h2 = ydim / NY
x, y = np.meshgrid(np.linspace(h1, xdim - h1, NX-1), np.linspace(h2, ydim - h2, NY-1))
F = f(x, y)
plt.imshow(F)
plt.show()

#solving the sparse system
flx = np.reshape(F, ((NX - 1) * (NY - 1), 1)) #lexicological ordering
u = np.linalg.solve(FDLaplacian(NX, NY), flx) #solving
u = np.reshape(u, ((NX - 1) , (NY - 1)))
plt.figure(2)
plt.clf()
plt.imshow(u)
plt.colorbar()
plt.show()



