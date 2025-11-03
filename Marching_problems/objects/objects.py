import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import identity

class Domain:
    def __init__(self, L, N, f, params_f, k, k_vectorized, params_k):
        self.grid = Grid(L, N)
        self.Forcing = Forcing(f, self.grid, params_f)
        self.F = self.Forcing.F.reshape((self.grid.Nx - 1) * (self.grid.Ny - 1))
        self.Kfunction = Kfunction(k_vectorized, self.grid, params_k)
        self.FV = FV(self.grid, k,  params_k)
        self.explicit_bounds = lambda self, x : (2 * x / self.grid.hx ** 2 + 2 * x / self.grid.hy ** 2)
        self.u = np.zeros_like(self.Forcing.F)
        self.integrator = Integrator(self.FV.A, self.u, self.F, 1)


    def _euler_call(self, T_0 = 0, T_1 = 1, N = 1):
        #assert stability
        self.dt = (T_1 - T_0) / N
        self.dt_max = self.explicit_bounds(self, np.max(self.Kfunction.K))
        if self.dt_max < self.dt:
            raise ValueError("unstable timestep")
        self.u = np.zeros_like(self.F)
        T = T_0
        self.integrator.reset()
        self.integrator.dt = self.dt
        print(self.integrator.dt)
        print(self.integrator.u)
        for i in range(N):
            self.integrator.explicit_euler()
            #T += self.dt
        self.u = self.integrator.u

    def _trapezoidal_call(self, T_0 = 0, T_1 = 1, N = 1):
        self.dt = (T_1 - T_0) / N
        #no stability requirements
        self.integrator.reset()
        self.integrator.dt = self.dt
        for i in range(N):
            self.integrator.implicit_trapezoid()
        self.u = self.integrator.u

    def get_Kfunction(self):
        plt.imshow(self.Kfunction.K)
        plt.show()

    def get_Forcing(self):
        print(self.Forcing.F)
        plt.imshow(self.Forcing.F)
        plt.show()

    def integrate(self, T_0 = 0, T_1 = 1, N = 1, method = 'trapezoidal'):
        if method == 'trapezoidal':
            self._trapezoidal_call(T_0, T_1, N)
        else:
            self._euler_call(T_0, T_1, N)


class WaveEquation:
    def __init__(self, L, N, f, params_f, k, k_vectorized, params_k, T_0 = 0, T_1 = 1):
        self.grid = Grid(L, N)
        self.Forcing = Forcing(f, self.grid, params_f)
        self.F = self.Forcing.F.reshape((self.grid.Nx - 1) * (self.grid.Ny - 1))
        self.Kfunction = Kfunction(k_vectorized, self.grid, params_k)
        self.FV = FV(self.grid, k, params_k)
        self.explicit_bounds = lambda self, x: (2 * x / self.grid.hx ** 2 + 2 * x / self.grid.hy ** 2)
        self._initial_conditions(T_0, T_1)
        self.T = 0
        self._update_central_scheme()
        self._update_central_scheme()

    def _initial_conditions(self, T_0 = 0, T_1 = 1):
        self.u = np.zeros_like(self.F)
        self.u_prime = np.zeros_like(self.F)
        # determine time step:
        self.dt_max = (2 / np.sqrt(self.explicit_bounds(self, np.max(self.Kfunction.K)))) * 0.5
        self.N = (T_1 - T_0) / self.dt_max
        self.dt = 0.01#self.dt_max
        #determine u1
        self.u1 = self.u + self.u_prime * self.dt + 0.5 * (-self.FV.A.dot(self.u) + self.F) * self.dt ** 2
        #determine A1
        self.A1 = 2 * identity(self.FV.A.shape[0], format=self.FV.A.format) - self.FV.A * self.dt ** 2
        print(f"Initialized Wave Equation scheme with dt = {self.dt}")

    def _update_central_scheme(self):
        buffer = self.u1.copy()
        #update F
        self.grid.t = self.T
        self.Forcing.update(self.grid)
        self.F = self.Forcing.F.reshape((self.grid.Nx - 1) * (self.grid.Ny - 1))
        self.u1 = -self.u + (self.A1.dot(self.u1) + self.F * self.dt ** 2)
        self.u = buffer
        self.T += self.dt
        print(f"Timestamp -> {self.T}")

    def step(self):
        self._update_central_scheme()

class Grid:
    def __init__(self, L, N):
       self.Lx, self.Ly = L
       self.Nx, self.Ny = N
       self.t = 0
       self.hx, self.hy = (self.Lx / self.Nx, self.Ly / self.Ny)
       self.i, self.j = np.meshgrid(np.arange(self.Nx - 1), np.arange(self.Ny - 1))


class Forcing:
    def __init__(self, f, grid, params):
        self.F = f(grid, *params)
        self.f = f
        self.params = params
    def update(self, grid):
        self.F = self.f(grid, *self.params)


class Kfunction:
    def __init__(self, k, grid, params):
        self.K = k(grid, params)


class FV:
    def __init__(self, grid, k, params):
        self.k = k
        self.grid = grid
        self.A = sc.sparse.dok_matrix(((self.grid.Nx - 1) * (self.grid.Ny - 1), (self.grid.Nx - 1) * (self.grid.Ny - 1)))
        self.get_global_matrix()

    def _coef_FV(self, grid, i, j):
        return (
            -self.k(grid, i - 0.5, j) / grid.hx ** 2,
            -self.k(grid, i, j - 0.5) / grid.hy ** 2,
             self.k(grid, i - 0.5, j) /grid.hx ** 2 + self.k(grid, i, j - 0.5) / grid.hy ** 2 + self.k(grid, i, j + 0.5) / grid.hy ** 2 + self.k(grid, i + 0.5, j) / grid.hx ** 2,
            -self.k(grid, i + 0.5, j) / grid.hx ** 2,
            -self.k(grid, i , j + 0.5) / grid.hy ** 2
        )

    def _get_coef_FV_latex(self, debug = True):
        self.Cof = np.empty((self.grid.Ny - 1, self.grid.Nx - 1), dtype=object)
        for j in range(1, self.grid.Ny):
          for i in range(1, self.grid.Nx):
              constants = self._coef_FV(self.grid, i, j)
              self.Cof[j - 1, i - 1] = constants
        self.Cof = self.Cof.reshape((self.grid.Ny - 1) *  (self.grid.Nx - 1))

    def _get_global_matrix_latex(self, debug = False):
        #fill the array
        row_len = (self.grid.Nx - 1) * (self.grid.Ny - 1)
        for j in range(row_len):
            row = j
            self.A[j, j] = self.Cof[row][2]
            i_u, j_u = j   - (self.grid.Nx - 1) * int(j /  (self.grid.Nx - 1)), int(j /  (self.grid.Nx - 1))
            if j > 0 :self.A[j, j - 1] = self.Cof[row][0]
            if j < (row_len- 1) : self.A[j, j + 1] = self.Cof[row][3]
            if j >= self.grid.Nx - 1 : self.A[j, j - self.grid.Nx + 1] = self.Cof[row][1]
            if j < (row_len - self.grid.Nx + 1) : self.A[j, j + self.grid.Nx - 1] = self.Cof[row][4]
            if(i_u == 0 and j > 0) : self.A[j, j - 1] = 0
            if(i_u == self.grid.Nx - 2 and j < (row_len - self.grid.Nx + 1)) : self.A[j, j + 1] = 0
            if(j_u == self.grid.Ny - 2 and j < (row_len - self.grid.Nx + 1)) : self.A[j, j +   self.grid.Nx - 1] = 0
            if(j_u == 0 ) : self.A[j, j -  self.grid.Nx + 1] = 0
        if self.A is None: return -1
        self.A = self.A.tocsr()
        return 1

    def get_global_matrix(self):
        self._get_coef_FV_latex()
        self._get_global_matrix_latex()
        return self.A

class Integrator:

    def __init__(self, A, u_0, f, dt):
        self.dt = dt
        self.A = A
        self.u = u_0
        self.F = f
        self.A1 = None
        self.A2 = None

    def explicit_euler(self):
        self.uprime = -self.A.dot(self.u) + self.F
        self.u = self.u + self.dt * self.uprime

    def implicit_trapezoid(self):
        if self.A1 is None:
            self.A1 = identity(self.A.shape[0], format = self.A.format) + self.A * self.dt / 2
            self.A2 = identity(self.A.shape[0], format = self.A.format) - self.A * self.dt / 2
        self.u = spsolve(self.A1, self.A2.dot(self.u) + self.F * self.dt)

    def reset(self):
        self.u = np.zeros_like(self.F)
        self.A1 =None
        self.A2 =None