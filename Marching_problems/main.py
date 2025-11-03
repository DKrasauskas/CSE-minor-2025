from Marching_problems.Problems.setup import *
from Marching_problems.Problems.waves.waves import *

w = Waves(L, N, launch_params_wave)
d = Diffusion(L, N, launch_params_diffusion)
d.get()

