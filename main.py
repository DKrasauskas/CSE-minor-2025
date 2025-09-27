import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve

NX =  4
NY = 4
grid = ((0, 10), (0, 5))



class Grid:
    def __init__(self, grid, NX, NY):
        self.NX = NX
        self.NY = NY
        self.grid = grid
        self.k = lambda i, j: 1 #+ 0.1 * (self.hx * i + self.hy * j + (self.hx * i) * (self.hy * j))
        self.k1 = lambda i, j: 1 + 0.1 * (self.hx * i + self.hy * j + (self.hx * i) * (self.hy * j))
        self._init_grid()
        self.K = self._coef_FV(self.i, self.j)
        self.A = sc.sparse.dok_matrix(((self.NX - 1) * (self.NY - 1), (self.NX - 1) * (self.NY - 1)))
        self.index = lambda i, j, row: (row , row)
        print(self.hx)
        print(self.hy)
        self._assemble_FD_global_matrix()

    def _init_grid(self):
        self.h = lambda maxV, minV, num : (maxV - minV) / num
        self.hx = self.h(self.grid[0][1], self.grid[0][0], self.NX)
        self.hy = self.h(self.grid[1][1], self.grid[1][0], self.NY)
        self.i, self.j = np.meshgrid(np.arange(self.NX - 1), np.arange(self.NY - 1))

    def _f(self, i, j, alpha=40):
        f_val = 0
        for k in range(1, 10):
            for z in range(1, 5):
                f_val += np.exp(-alpha * (i * self.hx - k) ** 2 - alpha * (j * self.hy - z) ** 2)
        return f_val

    def _assemble_FD_global_matrix(self, debug = False):
        self.DX = sc.sparse.dok_matrix((self.NX, self.NX - 1))
        self.DY = sc.sparse.dok_matrix((self.NY, self.NY - 1))
        for i in range(0, self.NX - 1):
            self.DX[i, i] = 1
            if i < (self.NX  -1): self.DX[i + 1, i] = -1
        for i in range(0, self.NY - 1):
            self.DY[i, i] = 1
            if i < (self.NY - 1): self.DY[i + 1, i] = -1
        self.Lxx = 1/self.hx ** 2 * (self.DX.T @ self.DX)
        self.Lyy = 1/self.hy ** 2 *(self.DY.T @ self.DY)
        print("Lyy =\n", self.Lyy.toarray())
        self.Ixx = sc.sparse.eye(self.NX - 1, self.NX - 1)
        self.Iyy = sc.sparse.eye(self.NY - 1, self.NY - 1)
        self.A1 = sc.sparse.kron(self.Iyy,  self.Lxx)
        self.A2 =  sc.sparse.kron(self.Lyy , self.Ixx)
        self.A_FD = (self.A1 + self.A2).tocsr()
        if debug:
            print("\n\n FD global matrix:")
            print(self.A_FD.toarray())

    def _coef_FV(self, i, j):
        return (
            -self.k(i - 0.5, j) / self.hx ** 2,
            -self.k(i, j - 0.5) / self.hy ** 2,
             self.k(i - 0.5, j) /self.hx ** 2 + self.k(i, j - 0.5) / self.hy ** 2 + self.k(i, j + 0.5) / self.hy ** 2 + self.k(i + 0.5, j) / self.hx ** 2,
            -self.k(i + 0.5, j) / self.hx ** 2,
            -self.k(i , j + 0.5) / self.hy ** 2
        )

    def get_coef_FV_latex(self, debug = False):
        row = 0
        self.Cof = np.empty((self.NY - 1, self.NX - 1), dtype=object)
        if debug: print("\n\nFV Coefficients")
        for j in range(1, self.NY):
          for i in range(1, self.NX):
              constants = self._coef_FV(i, j)
              self.Cof[j - 1, i - 1] = constants
              if debug: print(f"{constants[0]}u_{{{i - 1}, {j}}}  {constants[1]}u_{{{i}, {j - 1}}} +{constants[2]}u_{{{i}, {j}}}  {constants[3]}u_{{{i + 1}, {j}}}  {constants[4]}u_{{{i}, {j + 1}}}\\\\")
        self.Cof = self.Cof.reshape((self.NY - 1) *  (self.NX - 1))
        if not debug: return 0
        print("\n\nFV Coefficient Matrix")
        print(self.Cof)

    def get_global_matrix_latex(self, debug = False):
        #fill the array
        row_len = (self.NX - 1) * (self.NY - 1)
        for j in range((self.NX - 1) * (self.NY - 1)):
            row = j
            self.A[j, j] = self.Cof[row][2]
            i_u, j_u = j   - (self.NX - 1) * int(j /  (self.NX - 1)), int(j /  (self.NX - 1))
            if j > 0 :self.A[j, j - 1] = self.Cof[row][0]
            if j < (row_len- 1) : self.A[j, j + 1] = self.Cof[row][3]
            if j >= self.NX - 1 : self.A[j, j - self.NX + 1] = self.Cof[row][1]
            if j < (row_len - self.NX + 1) : self.A[j, j + self.NX - 1] = self.Cof[row][4]
            if(i_u == 0 and j > 0) : self.A[j, j - 1] = 0
            if(i_u == self.NX - 2 and j < (row_len - self.NX + 1)) : self.A[j, j + 1] = 0
            if(j_u == self.NY - 2 and j < (row_len - self.NX + 1)) : self.A[j, j +   self.NX - 1] = 0
            if(j_u == 0 ) : self.A[j, j -  self.NX + 1] = 0
        if self.A is None: return -1
        self.A = self.A.tocsr()
        if not debug: return 0
        self.B = self.A.toarray()
        print("\n\nglobal FV Matrix :")
        print(self.B)
        print("\n\nglobal FV Matrix Latex Style :")
        for i in range(0, len(self.B)):
            for j in  range(0, len(self.B[i])):
                print(self.B[i][j], end="&") if j != len(self.B) - 1 else print(self.B[i][j], end="\\\\")
            print()
        return 1

    def get_f_vector(self, debug = False):
        self.F = np.zeros((self.NX - 1) * (self.NY - 1))
        for j in range(1, self.NY):
            for i in range(1, self.NX):
                self.F[i -1 + (j - 1) * (self.NX -1)] = self._f(i, j)
        if not debug: return 0
        plt.imshow(self.F.reshape(self.NY-1, self.NX-1), extent=[grid[0][0], grid[0][1], grid[1][0], grid[1][1]])
        plt.colorbar()
        plt.show()

    def get_k_visuals(self):
        plt.imshow(np.vectorize(self.k)(self.i, self.j), extent=[grid[0][0], grid[0][1], grid[1][0], grid[1][1]])
        plt.colorbar()
        plt.show()


    def get_solution_FD(self, debug = False):
        self._assemble_FD_global_matrix(debug)
        self.get_f_vector(debug)
        self.u_FD = spsolve(self.A_FD, self.F).reshape(NY -1, NX-1)
        return self.u_FD

    def get_solution_FV(self, debug = False):
        self.get_coef_FV_latex(debug)
        self.get_global_matrix_latex(debug)
        self.get_f_vector(True)
        self.u_FV =  spsolve(self.A, self.F).reshape(NY -1, NX-1)
        return self.u_FV



G = Grid(grid, NX, NY)

G.get_k_visuals()
plt.imshow(G.get_solution_FD(debug=True), extent=[grid[0][0], grid[0][1], grid[1][0], grid[1][1]])
plt.colorbar()
plt.show()
plt.imshow(G.get_solution_FV(debug=True), extent=[grid[0][0], grid[0][1], grid[1][0], grid[1][1]])
plt.colorbar()
plt.show()

plt.imshow(G.A_FD.toarray())
plt.colorbar()
plt.show()



