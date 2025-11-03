import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from Marching_problems.objects.objects import WaveEquation, Domain
from Marching_problems.rendering.renderer import Renderer



class Waves:
    def __init__(self, L, N, params):
        self.L = L
        self.N = N
        self.params = params
    def get(self):
        self.wave = WaveEquation(*self.params)
        self.u =  self.wave.u.reshape((( self.wave.grid.Ny - 1), ( self.wave.grid.Nx - 1)))

        self.params = (
            self.wave.grid.hx / 2,
            self.wave.grid.Lx - self.wave.grid.hx / 2,
            self.wave.grid.Ly - self.wave.grid.hy / 2,
            self.wave.grid.hy / 2
        )
        r = Renderer(self.u, self.wave.Forcing.F, self.wave.Kfunction.K, self.params)
        for i in range(0, 100000):
            self.wave.step()
            if i % 2 == 0:
                self.u = self.wave.u.reshape(((self.wave.grid.Ny - 1), (self.wave.grid.Nx - 1)))
                # plt.clf()
                r.new_frame(self.u, self.wave.T, self.wave.Forcing.F)

        print("success")

class Diffusion:
    def __init__(self, L, N, params):
        self.L = L
        self.N = N
        self.params = params
    def get(self):
        self.domain = Domain(*self.params)
        print(f"computing Explicit for N  {1600}")
        self.domain.integrate(0, 1, 1600, 'euler')
        print("Done")
        self.u1 =  self.domain .u.copy()
        print(f"computing Explicit for N {3200}")
        self.domain.integrate(0, 1, 3200, 'euler')
        print("Done")
        self.u2 =  self.domain .u.copy()
        #generate implicit solutions:
        self.rms11 = []
        self.rms12 = []
        self.rms13 = []
        self.index = []
        for i in range(0, 8):
            self.N = pow(2, i)
            self.domain.integrate(0, 1, self.N, 'trapezoidal')
            self.u = self.domain.u
            #now compute RMS values:
            self.rms1 = np.sqrt(np.mean((self.u1 - self.u) ** 2))
            self.rms2 = np.sqrt(np.mean((self.u2 - self.u) ** 2))

            self.rms11.append(self.rms1)
            self.rms12.append(self.rms2)
            self.index.append(i)

            print(f"computing implicit for N = {self.N}")
        plt.semilogy(self.index, self.rms12, color = 'red', label = 'RMS u_2')
        plt.semilogy(self.index, self.rms11, color = 'blue', label = 'RMS u_1')
        plt.legend()
        plt.show()