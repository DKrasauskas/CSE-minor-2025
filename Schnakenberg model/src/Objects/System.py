import numpy as np

import scipy as sc
from scipy.sparse import coo_matrix, hstack, csr_matrix, bmat, dia
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import bicgstab
from scipy.sparse.linalg import eigs
import time

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from renderer import Renderer
from Objects.Laplacian import *

np.random.seed(0)
class System:
    """
        The System class implements the coupled PDE system and methods for solving it.


        Attributes:
            params (tuple): a tuple containing the constants of the problem = a, b, Du, Dv, LX, LY, NX, NY
            renderer (Object): a dependency injection of the renderer class for visualizat
    """
    def __init__(self, params, renderer):
        self.a, self.b, self.Du, self.Dv, self.LX, self.LY, self.NX, self.NY = params
        self.h = self.LX / self.NX
        self.A = Laplacian((self.NX, self.NY)).A * 1 / self.h ** 2
        #initialize initial conditions
        self.u = self.u0() * np.ones((self.NX, self.NY)) + 0.005 * (self.a * self.b ) * np.random.rand(self.NX, self.NY)
        self.v = self.v0() * np.ones((self.NX, self.NY))
        self.u = self.u.reshape((self.NX *self.NY), 1)
        self.v = self.v.reshape((self.NX *self.NY), 1)
        self.k = 6
        self.kI = self.k * sc.sparse.eye(*self.A.shape, format='csr')
        self.ka = self.a * self.k * sc.sparse.csr_matrix(np.ones((self.NX, self.NY)).reshape((self.NX *self.NY), 1))
        self.kb = self.b * self.k * sc.sparse.csr_matrix(np.ones((self.NX, self.NY)).reshape((self.NX *self.NY), 1))
        self.it_count = []
        self.it = []
        self.renderer = renderer
        self.I = None
        self.r = None
        self.wallclock = 0
        self.debug = False

    def u0(self,  a=0.1305, b=0.7695):
        """
        Initial conditions for the U variable
        """
        return a + b

    def v0(self, a=0.1305, b=0.7695):
        """
        Initial conditions for the V variable
        """
        return b / (a + b) ** 2

    def solve(self, integration_mode, debug = False):
        """
            For a given integration mode (And, if necessary, a custom timescale), solves the coupled system
            @params: integration_mode (string): the integration mode (EXPLICIT or IMPLICIT), can be passed as a boolean value
        """
        if(integration_mode == 1):
            self.explicit_integration()
        else:
            self.implicit_integration()

    def explicit_integration(self, T0 = 0, T1 = 20):
        plt.ion()
        self.dt = 0
        self._step()
        lambdas, vecs = self._asert_stability()
        self.dt =1.9 / abs(np.max(np.abs(lambdas)))
        self.T = T0
        print(f"\033[33mImplicit Integration\033[0m. Time step chosen -> \033[34m{self.dt}\033[0m.\n Maximum stable time step for stability -> \033[34m{2/ abs(np.max(np.abs(lambdas)))}\033[0m.")
        print(f"\033[31mIMPORTANT\033[0m: If the solver decides the timestep chosen was poor, it will auto-adjust for optimum convergence.")
        print(f"\033[31mIMPORTANT\033[0m: The maximum possible timestep should not be chosen; a dynamic change means that the solver will not converge / the timestep will have to be updated every itteration")
        print(f"which would be very expensive")
        for i in range(0, 100000):
            start = time.time()
            #passed arg has no effect
            self._step(0.01)
            end = time.time()
            self.T += self.dt
            if self.debug:
                self.wallclock += end - start
            #this is the dynamic tuning which can be enabled by uncommenting; In the problem, it was not asked to be implemented, so I leave it as is
            # if i % 400 == 0:
            #     lambdas, vecs = self._asert_stability()
            #     self.dt = 2 / abs(np.max(np.abs(lambdas)))
            #     print("new dt")
            #this controls the rendering speed
            if i % 200 == 0:
                self.renderer.dynamic_data1 = self.u.reshape((self.NX , self.NY))
                self.renderer.dynamic_data2 = self.v.reshape((self.NX , self.NY))
                self.renderer.dynamic_data3 = self.it_count
                self.renderer.it = self.it
                self.renderer.T = self.T
                self.renderer.dT = self.dt
                self.renderer.render()
                plt.draw()
                plt.pause(0.001)
            if self.T > T1:
                break

    def implicit_integration(self, T0 = 0, T1 = 20, render_speed = 4):
        plt.ion()
        self.dt = 0
        self._step()
        lambdas, vecs = self._asert_stability() # this generates the optimum timestep via Jacobian
        self.dt = .1
        self.T = 0
        print(f"\033[33mImplicit Integration\033[0m. Time step chosen ->  \033[34m{self.dt}\033[0m.\nMaximum stable time step for adequate convergence ->  \033[34m{10/ abs(np.max(np.abs(lambdas)))}\033[0m.")
        print(f"\033[31mIMPORTANT\033[0m: If the solver decides the timestep chosen was poor, it will auto-adjust for optimum convergence.")
        for i in range(0, 100000):
            #it is possible to adjust dt dynamically within the loop;
            start = time.time()
            self._NR_step()
            end = time.time()
            self.T += self.dt
            if self.debug:
                self.wallclock += end - start
            self.it.append(i)
            if i % render_speed == 0:
                self.renderer.dynamic_data1 = self.u.reshape((self.NX , self.NY))
                self.renderer.dynamic_data2 = self.v.reshape((self.NX , self.NY))
                self.renderer.T = self.T
                self.renderer.dT = self.dt
                self.renderer.it = self.it
                self.renderer.dynamic_data3 = self.it_count
                self.renderer.render()
                plt.draw()
                plt.pause(0.001)
            if self.T > T1:
                break

    def _step(self, dt = 0.1):
        """
        computes a single update of the variables (single step)
        """
        self.uv = self.u * self.v
        self.uu = self.u * self.u
        self.uuv = self.uu *  self.v
        self.fu = self.Du * self.A.dot(self.u) + self.ka + self.k * (-self.u + self.uuv)
        self.fv = self.Dv * self.A.dot(self.v) + self.kb + self.k * (-self.uuv)

        # compute growth rate:

        self.u += self.dt * self.fu
        self.v += self.dt * self.fv



    def _asert_stability(self):
        """
        computes the jacobian matrix and returns its eigenvalues / eigenvectors for stability estimate
        This method *can* be used for dynamically tuning the integration; However it is very expensive if applied at every step;
        Thus smarter techniques are necessary
        """
        # compute the jacobian matrix:
        self.diag1 = sc.sparse.diags((2 * self.k * self.uv).ravel(), format='csr')
        self.diag2 = sc.sparse.diags((self.k * self.uu).ravel(), format='csr')
        self.J = bmat(
            [[self.A * self.Du - self.kI + self.diag1, self.diag2], [-self.diag1, self.A * self.Dv - self.diag2]],
            format='csr'
        )

        return eigs(self.J, k=2, which='LM')

    def _get_jacobian(self):
        """
               computes the value of the J function (the definition can be found within the report)
        """
        # compute the jacobian matrix:
        self.diag1 = sc.sparse.diags((2 * self.k * self.uv).ravel(), format='csr')
        self.diag2 = sc.sparse.diags((self.k * self.uu).ravel(), format='csr')
        self.J = bmat(
            [[self.A * self.Du - self.kI + self.diag1, self.diag2], [-self.diag1, self.A * self.Dv - self.diag2]],
            format='csr'
        )

    def _get_g(self, dt=0.1):
        """
        computes the value of the G function (the definition can be found within the report)
        """
        self.uv = self.u * self.v
        self.uu = self.u * self.u
        self.uuv = self.uu * self.v
        self.fu = self.Du * self.A.dot(self.u) + self.ka + self.k * (-self.u + self.uuv)
        self.fv = self.Dv * self.A.dot(self.v) + self.kb + self.k * (-self.uuv)
        self.b = np.concatenate((
                self.u - self.u_init -self.dt * self.fu,
                self.v - self.v_init - self.dt * self.fv
        ))
    def _get_g_jacobian(self):
        """
        computes the Jacobian matrix of the G function (the definition can be found within the report)
        """
        # compute the jacobian matrix:
        self.diag1 = sc.sparse.diags((2 * self.k * self.uv).ravel(), format='csr')
        self.diag2 = sc.sparse.diags((self.k * self.uu).ravel(), format='csr')
        self.J = bmat(
            [[self.A * self.Du - self.kI + self.diag1, self.diag2], [-self.diag1, self.A * self.Dv - self.diag2]],
            format='csr'
        )
        if self.I is None:
            self.I = sc.sparse.eye(self.J.shape[0], self.J.shape[1], format='csr')
        return self.I - self.dt * self.J


    def _asert_stability_cheap(self):
        """
        A cheaper way to estimate the maximum eigenvalue; Not used in the problem, but could be if dynamic tuning becomes a problem
        """
        fu_u_max = np.max(np.abs(-self.k + 2 * self.k * self.u * self.v))  # bound ∂R_u/∂u
        fu_v_max = np.max(np.abs(self.k * self.u ** 2))  # bound ∂R_u/∂v
        fv_u_max = np.max(np.abs(-2 * self.k * self.u * self.v))  # bound ∂R_v/∂u
        fv_v_max = np.max(np.abs(-self.k * self.u ** 2))
        return fu_u_max + fu_v_max+ fv_u_max+ fv_v_max

    def _NR_step(self):
        """
        The coolest of methods -> utilizes the NR method to perform an implicit timestep
        The method is tuned dynamically; if the NR step struggles to find a solution within a reasonable
        amount of steps, the correction clause is invoked which generates a new, optimum time-step via self.assert_stability()
        :return:
        """
        self.u_init = self.u.copy()
        self.v_init = self.v.copy()
        correction = False
        for i in range(0, 100):
            self._get_g()
            jacobian = self._get_g_jacobian()
            if self.r is not None:
                self.r, info = bicgstab(jacobian, self.b, x0=self.r, maxiter=100) #this may not be always optimum; for adequate choices of thimestep, faster than spsolve
            else:
                self.r = spsolve(jacobian, self.b)
            deltau = self.r[:len(self.u)]
            deltau = deltau.reshape(-1, 1)
            residual = np.sqrt(np.dot(self.r, self.r))
            deltav = self.r[len(self.u):]
            deltav = deltav.reshape(-1, 1)
            #once solved, update:
            self.u -= deltau
            self.v -= deltav
            if residual < 1e-3:
                self.it_count.append(i)
                break
            if i > 3 and not correction:
                self.u = self.u_init.copy()
                self.v = self.v_init.copy()
                lambdas, vec = self._asert_stability() #update the timestep to an optimum / stable one
                self.dt = 10 / np.max(np.abs(lambdas)) #factor of 10 chosen from experimentation for best performance
                correction = True
                print(f"\033[31mCondition Number High\033[0m: Adjusting dt to  {self.dt}")
                print(f"\033[31mWarning\033[0m: Runtime is expected to increase. For faster simulation, make a better dt choice {self.dt}")
                print(f"\033[31mWarning\033[0m: Solver cannot risk stability in choosing a higher dt itself. {self.dt}")
                i = 0






